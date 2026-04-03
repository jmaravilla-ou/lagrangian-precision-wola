--------------------------------------------------------------------------------
-- WOLA_Filterbank.vhd                    (TCAS-II submission - full filterbank)
--
-- Top-level entity.  Wires the five components of the variable-precision
-- WOLA filterbank described in Section IV of the paper:
--
--   ┌──────────┐  sub_re/im/mag   ┌───────────────────┐
--   │  WOLA    │ ───────────────▶  │  Precision        │
--   │ Analysis │   frame_valid     │  Controller (FSM) │
--   └──────────┘                   └────────┬──────────┘
--                                     mask  │  band_index
--                                  ┌────────▼──────────┐
--        sub_re/im ──────────────▶ │ Gain Application  │ ◄── Gain LUT
--                                  │ (Isolation + MAC) │     (per-band)
--                                  └────────┬──────────┘
--                                     proc_re/im
--                                  ┌────────▼──────────┐
--                                  │  WOLA             │
--                                  │  Synthesis        │ ──▶ pcm_out
--                                  └───────────────────┘
--
-- Gain-application datapath
-- ─────────────────────────
-- For each subband k the Precision Controller drives:
--   • mask_out : 24-bit AND mask with upper B_k bits set
--   • band_index : selects which band's data to process
--
-- The MAC computes  gain_k × subband_re_k  and  gain_k × subband_im_k
-- (two passes per band, serialised) with operand-isolation applied to
-- both operands.  The processed re/im are then forwarded to synthesis.
--
-- The per-band gain vector defaults to unity (0x400000 in Q1.22) and can
-- be overwritten through the gain_we / gain_addr / gain_data port.
-- In a full hearing-aid SoC this gain vector would be updated by the
-- fitting algorithm running on a companion processor.
--
-- Target: Xilinx Artix-7 XC7A200T, 100 MHz, Vivado 2025.2 WebPACK
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity WOLA_Filterbank is
    port (
        clk           : in  std_logic;
        rst           : in  std_logic;
        -- Audio input
        pcm_in        : in  std_logic_vector(23 downto 0);
        pcm_in_valid  : in  std_logic;
        -- Precision allocation (loaded at startup)
        alloc_we      : in  std_logic;
        alloc_addr    : in  std_logic_vector(5 downto 0);
        alloc_data    : in  std_logic_vector(4 downto 0);
        -- Per-band gain vector (loaded at startup, Q1.22 signed)
        gain_we       : in  std_logic;
        gain_addr     : in  std_logic_vector(5 downto 0);
        gain_data     : in  std_logic_vector(23 downto 0);
        -- Audio output
        pcm_out       : out std_logic_vector(23 downto 0);
        pcm_out_valid : out std_logic
    );
end WOLA_Filterbank;

architecture Structural of WOLA_Filterbank is

    constant W  : integer := 24;
    constant NB : integer := 64;

    ---------------------------------------------------------------------------
    -- Component declarations
    ---------------------------------------------------------------------------
    component WOLA_Analysis is
        generic (WIDTH:integer:=24; FRAME_LEN:integer:=128;
                 HOP_SIZE:integer:=64; NUM_BANDS:integer:=64);
        port (clk,rst:in std_logic;
              pcm_in:in std_logic_vector(WIDTH-1 downto 0);
              pcm_in_valid:in std_logic;
              sub_re,sub_im,sub_mag:out std_logic_vector(WIDTH-1 downto 0);
              sub_index:out std_logic_vector(5 downto 0);
              sub_valid,frame_valid:out std_logic;
              frame_complete:out std_logic);
    end component;

    component Precision_Controller is
        generic (WIDTH:integer:=24; NUM_BANDS:integer:=64);
        port (clk,rst:in std_logic;
              alloc_we:in std_logic;
              alloc_addr:in std_logic_vector(5 downto 0);
              alloc_data:in std_logic_vector(4 downto 0);
              frame_valid:in std_logic;
              frame_done:out std_logic;
              band_index:out std_logic_vector(5 downto 0);
              mask_out:out std_logic_vector(WIDTH-1 downto 0);
              mac_en,band_valid:out std_logic;
              band_ack:in std_logic);
    end component;

    component Operand_Isolation is
        generic (WIDTH:integer:=24);
        port (a_raw,b_raw,mask:in std_logic_vector(WIDTH-1 downto 0);
              a_masked,b_masked:out std_logic_vector(WIDTH-1 downto 0));
    end component;

    component Precision_MAC is
        generic (WIDTH:integer:=24);
        port (clk,rst,en:in std_logic;
              a_in,b_in:in std_logic_vector(WIDTH-1 downto 0);
              res_out:out std_logic_vector(2*WIDTH-1 downto 0));
    end component;

    component WOLA_Synthesis is
        generic (WIDTH:integer:=24; FRAME_LEN:integer:=128;
                 HOP_SIZE:integer:=64; NUM_BANDS:integer:=64);
        port (clk,rst:in std_logic;
              sub_re,sub_im:in std_logic_vector(WIDTH-1 downto 0);
              sub_index:in std_logic_vector(5 downto 0);
              sub_valid,frame_start:in std_logic;
              pcm_out:out std_logic_vector(WIDTH-1 downto 0);
              pcm_out_valid:out std_logic);
    end component;

    ---------------------------------------------------------------------------
    -- Analysis outputs
    ---------------------------------------------------------------------------
    signal a_re, a_im, a_mag : std_logic_vector(W-1 downto 0);
    signal a_idx              : std_logic_vector(5 downto 0);
    signal a_val, a_fval      : std_logic;
    signal a_fcomplete        : std_logic;

    ---------------------------------------------------------------------------
    -- Subband storage (registered from analysis)
    ---------------------------------------------------------------------------
    type sub_arr_t is array (0 to NB-1) of std_logic_vector(W-1 downto 0);
    signal sr_mem, si_mem : sub_arr_t := (others => (others => '0'));

    ---------------------------------------------------------------------------
    -- Per-band gain LUT  (Q1.22 signed - unity = 0x400000)
    ---------------------------------------------------------------------------
    type gain_arr_t is array (0 to NB-1) of std_logic_vector(W-1 downto 0);
    signal gain_lut : gain_arr_t := (others => x"400000"); -- default unity

    ---------------------------------------------------------------------------
    -- Precision Controller
    ---------------------------------------------------------------------------
    signal pc_done                   : std_logic;
    signal pc_idx                    : std_logic_vector(5 downto 0);
    signal pc_mask                   : std_logic_vector(W-1 downto 0);
    signal pc_en, pc_bval            : std_logic;
    signal pc_band_ack               : std_logic := '0';

    ---------------------------------------------------------------------------
    -- Gain-application sub-FSM
    --
    -- For each band the Precision Controller asserts band_valid for one clock.
    -- We need TWO multiply passes:  gain × re  then  gain × im.
    -- A small sub-FSM (G_IDLE → G_RE → G_CAP_RE → G_IM → G_CAP_IM → G_STORE)
    -- controls the mux into the Operand_Isolation + MAC pipeline.
    --
    -- Handshake: The PC enters S_WAIT after S_EXECUTE and holds until it
    -- sees pc_band_ack = '1', which the sub-FSM asserts in G_STORE.
    -- This ensures the PC never issues a new band_valid while the gain
    -- pipeline is still processing the previous band.
    ---------------------------------------------------------------------------
    type gst_t is (G_IDLE, G_RE, G_CAP_RE, G_IM, G_CAP_IM, G_STORE);
    signal gst : gst_t := G_IDLE;

    signal band_latch   : unsigned(5 downto 0)       := (others => '0');
    signal mask_latch   : std_logic_vector(W-1 downto 0) := (others => '0');
    signal gain_latch   : std_logic_vector(W-1 downto 0) := (others => '0');

    signal iso_a, iso_b : std_logic_vector(W-1 downto 0) := (others => '0');
    signal iso_am, iso_bm : std_logic_vector(W-1 downto 0) := (others => '0');
    signal mac_en_i     : std_logic := '0';
    signal mac_res      : std_logic_vector(2*W-1 downto 0);

    -- Processed re/im per band (Q1.22 × Q0.23 → take upper 24 of 48)
    type proc_arr_t is array (0 to NB-1) of std_logic_vector(W-1 downto 0);
    signal proc_re_mem, proc_im_mem : proc_arr_t := (others => (others => '0'));
    signal proc_re_cap : std_logic_vector(W-1 downto 0) := (others => '0');

    ---------------------------------------------------------------------------
    -- Synthesis feed sequencer
    ---------------------------------------------------------------------------
    signal sf_active  : std_logic := '0';
    signal sf_idx     : unsigned(5 downto 0) := (others => '0');
    signal sf_val     : std_logic := '0';
    signal sf_fstart  : std_logic := '0';

    signal s_re, s_im : std_logic_vector(W-1 downto 0) := (others => '0');
    signal s_idx      : std_logic_vector(5 downto 0) := (others => '0');

begin

    ---------------------------------------------------------------------------
    -- 1.  WOLA Analysis
    ---------------------------------------------------------------------------
    U_ANA : WOLA_Analysis
        generic map (WIDTH=>W, FRAME_LEN=>128, HOP_SIZE=>64, NUM_BANDS=>NB)
        port map (clk=>clk, rst=>rst,
                  pcm_in=>pcm_in, pcm_in_valid=>pcm_in_valid,
                  sub_re=>a_re, sub_im=>a_im, sub_mag=>a_mag,
                  sub_index=>a_idx, sub_valid=>a_val, frame_valid=>a_fval,
                  frame_complete=>a_fcomplete);

    -- Buffer analysis output into per-band RAMs --------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if a_val = '1' then
                sr_mem(to_integer(unsigned(a_idx))) <= a_re;
                si_mem(to_integer(unsigned(a_idx))) <= a_im;
            end if;
        end if;
    end process;

    ---------------------------------------------------------------------------
    -- 2.  Precision Controller
    ---------------------------------------------------------------------------
    U_PC : Precision_Controller
        generic map (WIDTH=>W, NUM_BANDS=>NB)
        port map (clk=>clk, rst=>rst,
                  alloc_we=>alloc_we, alloc_addr=>alloc_addr,
                  alloc_data=>alloc_data,
                  frame_valid=>a_fcomplete, frame_done=>pc_done,
                  band_index=>pc_idx, mask_out=>pc_mask,
                  mac_en=>pc_en, band_valid=>pc_bval,
                  band_ack=>pc_band_ack);

    ---------------------------------------------------------------------------
    -- Per-band gain LUT write
    ---------------------------------------------------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if gain_we = '1' then
                gain_lut(to_integer(unsigned(gain_addr))) <= gain_data;
            end if;
        end if;
    end process;

    ---------------------------------------------------------------------------
    -- 3 + 4.  Gain-application sub-FSM  (Operand Isolation + MAC)
    --
    -- Each band takes 5 clocks: G_RE → G_CAP_RE → G_IM → G_CAP_IM → G_STORE.
    -- The Precision Controller waits in S_WAIT for pc_band_ack before
    -- advancing to the next band, so there is no timing conflict.
    --
    -- Timing:  G_RE drives operands → MAC captures product → G_CAP_RE reads
    -- mac_res (1-cycle MAC latency).  Same for G_IM / G_CAP_IM.
    ---------------------------------------------------------------------------
    process(clk)
        variable ki : integer;
    begin
        if rising_edge(clk) then
            if rst = '1' then
                gst <= G_IDLE; mac_en_i <= '0'; pc_band_ack <= '0';
            else
                mac_en_i    <= '0';
                pc_band_ack <= '0';

                case gst is
                when G_IDLE =>
                    if pc_bval = '1' then
                        band_latch <= unsigned(pc_idx);
                        mask_latch <= pc_mask;
                        gain_latch <= gain_lut(to_integer(unsigned(pc_idx)));
                        gst <= G_RE;
                    end if;

                when G_RE =>
                    -- Drive Operand A = gain_k, Operand B = sub_re_k
                    ki := to_integer(band_latch);
                    iso_a    <= gain_latch;
                    iso_b    <= sr_mem(ki);
                    mac_en_i <= '1';
                    gst      <= G_CAP_RE;

                when G_CAP_RE =>
                    -- MAC result available this cycle (1-cycle latency)
                    -- Product is Q1.22 × Q0.23 = Q1.45 (48 bits).
                    -- Take bits [45:22] as Q1.22 result.
                    proc_re_cap <= mac_res(2*W-3 downto W-2);
                    -- Drive IM operands
                    ki := to_integer(band_latch);
                    iso_a    <= gain_latch;
                    iso_b    <= si_mem(ki);
                    mac_en_i <= '1';
                    gst      <= G_IM;

                when G_IM =>
                    -- Store RE result
                    proc_re_mem(to_integer(band_latch)) <= proc_re_cap;
                    gst <= G_CAP_IM;

                when G_CAP_IM =>
                    -- Capture IM product
                    proc_im_mem(to_integer(band_latch)) <= mac_res(2*W-3 downto W-2);
                    gst <= G_STORE;

                when G_STORE =>
                    pc_band_ack <= '1';   -- signal PC to advance to next band
                    gst <= G_IDLE;

                end case;
            end if;
        end if;
    end process;

    ---------------------------------------------------------------------------
    -- Operand Isolation instance
    ---------------------------------------------------------------------------
    U_ISO : Operand_Isolation
        generic map (WIDTH=>W)
        port map (a_raw=>iso_a, b_raw=>iso_b, mask=>mask_latch,
                  a_masked=>iso_am, b_masked=>iso_bm);

    ---------------------------------------------------------------------------
    -- Precision MAC instance
    ---------------------------------------------------------------------------
    U_MAC : Precision_MAC
        generic map (WIDTH=>W)
        port map (clk=>clk, rst=>rst, en=>mac_en_i,
                  a_in=>iso_am, b_in=>iso_bm, res_out=>mac_res);

    ---------------------------------------------------------------------------
    -- Synthesis feed sequencer
    --   After the Precision Controller signals frame_done (all 64 bands
    --   processed through the gain stage), stream the processed re/im
    --   to the WOLA_Synthesis block at one band per clock.
    ---------------------------------------------------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if rst = '1' then
                sf_active <= '0'; sf_idx <= (others=>'0');
                sf_val <= '0'; sf_fstart <= '0';
            else
                sf_val    <= '0';
                sf_fstart <= '0';

                if pc_done = '1' and sf_active = '0' then
                    sf_active <= '1';
                    sf_idx    <= (others=>'0');
                    sf_fstart <= '1';
                end if;

                if sf_active = '1' then
                    s_re  <= proc_re_mem(to_integer(sf_idx));
                    s_im  <= proc_im_mem(to_integer(sf_idx));
                    s_idx <= std_logic_vector(sf_idx);
                    sf_val <= '1';

                    if sf_idx = to_unsigned(NB-1, 6) then
                        sf_active <= '0';
                    else
                        sf_idx <= sf_idx + 1;
                    end if;
                end if;
            end if;
        end if;
    end process;

    ---------------------------------------------------------------------------
    -- 5.  WOLA Synthesis
    ---------------------------------------------------------------------------
    U_SYN : WOLA_Synthesis
        generic map (WIDTH=>W, FRAME_LEN=>128, HOP_SIZE=>64, NUM_BANDS=>NB)
        port map (clk=>clk, rst=>rst,
                  sub_re=>s_re, sub_im=>s_im, sub_index=>s_idx,
                  sub_valid=>sf_val, frame_start=>sf_fstart,
                  pcm_out=>pcm_out, pcm_out_valid=>pcm_out_valid);

end Structural;