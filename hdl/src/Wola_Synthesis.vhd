--------------------------------------------------------------------------------
-- WOLA_Synthesis_v2.vhd                  (TCAS-II submission - full filterbank)
--
-- Synthesisable WOLA synthesis back-end using Xilinx xfft_0 IP (v9.1.15)
-- configured for inverse FFT.
-- Replaces the behavioural real-arithmetic IFFT in Wola_Synthesis.vhd.
--
-- PIPELINE:
--   Collect 64 bins → mirror to 128-pt conjugate-symmetric spectrum
--   → xfft_0 (IFFT, 128pt) → synthesis Hann window → overlap-add → PCM out
--
-- IMPORTANT - TWO XFFT_0 INSTANCES:
--   The analysis FFT and synthesis IFFT each need their own xfft_0 instance
--   because they run concurrently and have independent state.
--   This file instantiates a SECOND xfft_0, named xfft_1_inst, configured
--   as IFFT via the config channel (FWD_INV bit = 1).
--
--   Add xfft_0.vhd to your project ONCE; it is instantiated twice:
--     WOLA_Analysis_v2.vhd  →  U_FFT  (forward, FWD_INV=0)
--     WOLA_Synthesis_v2.vhd →  U_IFFT (inverse, FWD_INV=1)
--
-- SCALING:
--   IFFT uses SCALE_SCH = 8'b00000000 (no scaling) to complement the
--   analysis FFT's divide-by-8. The WOLA round-trip gain is:
--     analysis_scale(1/8) × IFFT_scale(1) × OLA_normalisation ≈ 1
--   The Hann window overlap-add with 50% overlap normalises the gain
--   to unity for the baseband signal.
--
-- CONFIG WORD for IFFT:
--   bits[8:1] = SCALE_SCH = 0x00 (no scaling)
--   bit[0]    = FWD_INV   = 1    (inverse)
--   Full 16-bit = 0x0001
--
-- OUTPUT:
--   Emits HOP_SIZE (64) PCM samples per frame after overlap-add.
--   pcm_out_valid pulses once per output sample.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use IEEE.MATH_REAL.ALL;

entity WOLA_Synthesis is
    generic (
        WIDTH      : integer := 24;
        FRAME_LEN  : integer := 128;
        HOP_SIZE   : integer := 64;
        NUM_BANDS  : integer := 64
    );
    port (
        clk           : in  std_logic;
        rst           : in  std_logic;
        sub_re        : in  std_logic_vector(WIDTH-1 downto 0);
        sub_im        : in  std_logic_vector(WIDTH-1 downto 0);
        sub_index     : in  std_logic_vector(5 downto 0);
        sub_valid     : in  std_logic;
        frame_start   : in  std_logic;
        pcm_out       : out std_logic_vector(WIDTH-1 downto 0);
        pcm_out_valid : out std_logic
    );
end WOLA_Synthesis;

architecture RTL of WOLA_Synthesis is

    ---------------------------------------------------------------------------
    -- xfft_0 component declaration (same IP, configured as IFFT)
    ---------------------------------------------------------------------------
    component xfft_0 is
        port (
            aclk                        : in  std_logic;
            s_axis_config_tdata         : in  std_logic_vector(15 downto 0);
            s_axis_config_tvalid        : in  std_logic;
            s_axis_config_tready        : out std_logic;
            s_axis_data_tdata           : in  std_logic_vector(47 downto 0);
            s_axis_data_tvalid          : in  std_logic;
            s_axis_data_tready          : out std_logic;
            s_axis_data_tlast           : in  std_logic;
            m_axis_data_tdata           : out std_logic_vector(47 downto 0);
            m_axis_data_tvalid          : out std_logic;
            m_axis_data_tready          : in  std_logic;
            m_axis_data_tlast           : out std_logic;
            event_frame_started         : out std_logic;
            event_tlast_unexpected      : out std_logic;
            event_tlast_missing         : out std_logic;
            event_status_channel_halt   : out std_logic;
            event_data_in_channel_halt  : out std_logic;
            event_data_out_channel_halt : out std_logic
        );
    end component;

    ---------------------------------------------------------------------------
    -- Synthesis Hann ROM (unsigned, same coefficients as analysis)
    ---------------------------------------------------------------------------
    type rom_t is array (0 to FRAME_LEN-1) of unsigned(WIDTH-1 downto 0);
    function init_hann return rom_t is
        variable r : rom_t;
    begin
        for n in 0 to FRAME_LEN-1 loop
            r(n) := to_unsigned(integer(round(
                0.5*(1.0 - cos(2.0*MATH_PI*real(n)/real(FRAME_LEN)))
                * real(2**(WIDTH-1) - 1))), WIDTH);
        end loop;
        return r;
    end function;
    constant HANN_ROM : rom_t := init_hann;

    ---------------------------------------------------------------------------
    -- Spectrum collection buffers (one-sided input, full 128-bin for IFFT)
    ---------------------------------------------------------------------------
    type cbuf_t is array (0 to FRAME_LEN-1) of signed(WIDTH-1 downto 0);
    signal spec_re : cbuf_t := (others => (others => '0'));
    signal spec_im : cbuf_t := (others => (others => '0'));

    signal bcnt         : unsigned(5 downto 0) := (others => '0');
    signal collect_done : std_logic := '0';

    ---------------------------------------------------------------------------
    -- FSM states
    -- COLLECT  : accept 64 one-sided bins from upstream
    -- MIRROR   : build conjugate-symmetric 128-point spectrum (64 cycles)
    -- FEED_CFG : send IFFT config to xfft_0
    -- FEED     : stream 128 bins to xfft_0
    -- IFFT_WAIT: wait for IFFT output (all 128 samples)
    -- WIN_OLA  : apply synthesis window + overlap-add (FRAME_LEN cycles)
    -- OUTPUT   : emit HOP_SIZE samples
    ---------------------------------------------------------------------------
    type state_t is (S_COLLECT, S_MIRROR, S_FEED_CFG, S_FEED,
                     S_IFFT_WAIT, S_WIN_OLA, S_OUTPUT);
    signal state : state_t := S_COLLECT;

    signal mirror_idx : unsigned(5 downto 0) := (others => '0');
    signal feed_idx   : unsigned(6 downto 0) := (others => '0');

    -- IFFT input signals
    signal ifft_din_tdata  : std_logic_vector(47 downto 0) := (others => '0');
    signal ifft_din_tvalid : std_logic := '0';
    signal ifft_din_tready : std_logic;
    signal ifft_din_tlast  : std_logic := '0';

    -- IFFT config
    signal ifft_cfg_tdata  : std_logic_vector(15 downto 0) := (others => '0');
    signal ifft_cfg_tvalid : std_logic := '0';
    signal ifft_cfg_tready : std_logic;

    -- IFFT output
    signal ifft_dout_tdata  : std_logic_vector(47 downto 0);
    signal ifft_dout_tvalid : std_logic;
    signal ifft_dout_tlast  : std_logic;

    -- IFFT time-domain output buffer
    signal ifft_buf  : cbuf_t := (others => (others => '0'));
    signal ifft_idx  : unsigned(6 downto 0) := (others => '0');
    signal ifft_done : std_logic := '0';

    -- Windowed current frame and OLA tail
    type ola_t is array (0 to FRAME_LEN-1) of signed(WIDTH downto 0);
    signal ola_curr : ola_t := (others => (others => '0'));
    signal ola_prev : ola_t := (others => (others => '0'));

    -- Window + OLA state
    signal win_idx    : unsigned(6 downto 0) := (others => '0');

    -- Output sequencer
    signal out_idx    : unsigned(6 downto 0) := (others => '0');
    signal oval       : std_logic := '0';

begin

    ---------------------------------------------------------------------------
    -- xfft_0 instantiation - configured as IFFT
    -- Config: SCALE_SCH=0x00 (no scaling), FWD_INV=1 (inverse)
    -- 16-bit config word = 0x0001
    ---------------------------------------------------------------------------
    ifft_cfg_tdata <= x"0001";   -- IFFT, no scaling

    U_IFFT : xfft_0
        port map (
            aclk                        => clk,
            s_axis_config_tdata         => ifft_cfg_tdata,
            s_axis_config_tvalid        => ifft_cfg_tvalid,
            s_axis_config_tready        => ifft_cfg_tready,
            s_axis_data_tdata           => ifft_din_tdata,
            s_axis_data_tvalid          => ifft_din_tvalid,
            s_axis_data_tready          => ifft_din_tready,
            s_axis_data_tlast           => ifft_din_tlast,
            m_axis_data_tdata           => ifft_dout_tdata,
            m_axis_data_tvalid          => ifft_dout_tvalid,
            m_axis_data_tready          => '1',
            m_axis_data_tlast           => ifft_dout_tlast,
            event_frame_started         => open,
            event_tlast_unexpected      => open,
            event_tlast_missing         => open,
            event_status_channel_halt   => open,
            event_data_in_channel_halt  => open,
            event_data_out_channel_halt => open
        );

    ---------------------------------------------------------------------------
    -- Main FSM
    ---------------------------------------------------------------------------
    p_main: process(clk)
        variable ki   : integer;
        variable wp   : signed(2*(WIDTH+1)-1 downto 0);
        variable ws   : signed(WIDTH-1 downto 0);
        variable osum : signed(WIDTH downto 0);
        variable re_v, im_v : signed(WIDTH-1 downto 0);
    begin
        if rising_edge(clk) then
            if rst = '1' then
                state          <= S_COLLECT;
                bcnt           <= (others => '0');
                collect_done   <= '0';
                mirror_idx     <= (others => '0');
                feed_idx       <= (others => '0');
                ifft_idx       <= (others => '0');
                win_idx        <= (others => '0');
                out_idx        <= (others => '0');
                ifft_done      <= '0';
                oval           <= '0';
                ifft_din_tvalid<= '0';
                ifft_din_tlast <= '0';
                ifft_cfg_tvalid<= '0';
            else
                collect_done   <= '0';
                ifft_done      <= '0';
                oval           <= '0';

                case state is

                -- ------------------------------------------------------------
                -- S_COLLECT: accept one-sided bins 0..63 from gain stage
                -- Also mirror bin 0 (DC, real-valued) on arrival
                -- ------------------------------------------------------------
                when S_COLLECT =>
                    ifft_din_tvalid <= '0';

                    if frame_start = '1' then
                        bcnt <= (others => '0');
                        -- Zero the Nyquist bin (index 64)
                        spec_re(NUM_BANDS) <= (others => '0');
                        spec_im(NUM_BANDS) <= (others => '0');
                    end if;

                    if sub_valid = '1' then
                        ki := to_integer(unsigned(sub_index));
                        spec_re(ki) <= signed(sub_re);
                        spec_im(ki) <= signed(sub_im);
                        bcnt <= bcnt + 1;

                        if bcnt = to_unsigned(NUM_BANDS - 1, 6) then
                            collect_done <= '1';
                            state        <= S_MIRROR;
                            mirror_idx   <= to_unsigned(1, 6);
                        end if;
                    end if;

                -- ------------------------------------------------------------
                -- S_MIRROR: fill bins 65..127 with conjugate of bins 63..1
                -- X[N-k] = conj(X[k]) for k = 1..63
                -- Takes 63 cycles
                -- ------------------------------------------------------------
                when S_MIRROR =>
                    ki := to_integer(mirror_idx);
                    spec_re(FRAME_LEN - ki) <=  spec_re(ki);
                    spec_im(FRAME_LEN - ki) <= -spec_im(ki);

                    if mirror_idx = to_unsigned(NUM_BANDS - 1, 6) then
                        -- Done mirroring; send config to IFFT
                        ifft_cfg_tvalid <= '1';
                        state           <= S_FEED_CFG;
                    else
                        mirror_idx <= mirror_idx + 1;
                    end if;

                -- ------------------------------------------------------------
                -- S_FEED_CFG: wait for IFFT to accept config
                -- ------------------------------------------------------------
                when S_FEED_CFG =>
                    if ifft_cfg_tready = '1' then
                        ifft_cfg_tvalid <= '0';
                        feed_idx        <= (others => '0');
                        state           <= S_FEED;
                    end if;

                -- ------------------------------------------------------------
                -- S_FEED: stream 128 complex bins to xfft_0 (IFFT mode)
                -- Data packing: [47:24]=Re, [23:0]=Im
                -- Respect tready backpressure
                -- ------------------------------------------------------------
                when S_FEED =>
                    if ifft_din_tready = '1' or ifft_din_tvalid = '0' then
                        ki := to_integer(feed_idx);
                        -- Pack Re into [47:24], Im into [23:0]
                        ifft_din_tdata  <= std_logic_vector(spec_re(ki))
                                         & std_logic_vector(spec_im(ki));
                        ifft_din_tvalid <= '1';
                        ifft_din_tlast  <= '0';

                        if feed_idx = to_unsigned(FRAME_LEN - 1, 7) then
                            ifft_din_tlast <= '1';
                            ifft_idx       <= (others => '0');
                            state          <= S_IFFT_WAIT;
                        else
                            feed_idx <= feed_idx + 1;
                        end if;
                    end if;

                -- ------------------------------------------------------------
                -- S_IFFT_WAIT: collect 128 IFFT output samples
                -- Only Re part is meaningful (imaginary is near-zero for
                -- conjugate-symmetric input). Store Re to ifft_buf.
                -- ------------------------------------------------------------
                when S_IFFT_WAIT =>
                    ifft_din_tvalid <= '0';
                    ifft_din_tlast  <= '0';

                    if ifft_dout_tvalid = '1' then
                        -- Re is in [47:24]
                        re_v := signed(ifft_dout_tdata(47 downto 24));
                        ifft_buf(to_integer(ifft_idx)) <= re_v;

                        if ifft_idx = to_unsigned(FRAME_LEN - 1, 7) then
                            ifft_done <= '1';
                            win_idx   <= (others => '0');
                            state     <= S_WIN_OLA;
                        else
                            ifft_idx <= ifft_idx + 1;
                        end if;
                    end if;

                -- ------------------------------------------------------------
                -- S_WIN_OLA: apply synthesis Hann window to all 128 samples
                -- + compute OLA sum (current + prev tail)
                -- Takes FRAME_LEN cycles
                -- Product: Q1.22 × Q0.23 (unsigned Hann) → take [2W-2:W-1]
                -- OLA buffer is WIDTH+1 bits to hold sum without overflow
                -- ------------------------------------------------------------
                when S_WIN_OLA =>
                    wp := resize(ifft_buf(to_integer(win_idx)), WIDTH+1) *
                          signed('0' & std_logic_vector(HANN_ROM(to_integer(win_idx))));
                    ws := wp(2*WIDTH-1 downto WIDTH);
                    ola_curr(to_integer(win_idx)) <= resize(ws, WIDTH+1);

                    if win_idx = to_unsigned(FRAME_LEN - 1, 7) then
                        out_idx <= (others => '0');
                        state   <= S_OUTPUT;
                    else
                        win_idx <= win_idx + 1;
                    end if;

                -- ------------------------------------------------------------
                -- S_OUTPUT: emit HOP_SIZE overlap-added samples
                -- PCM[n] = ola_curr[n] + ola_prev[n + HOP_SIZE]
                -- Saturate to WIDTH bits
                -- After emitting, save ola_curr as ola_prev for next frame
                -- ------------------------------------------------------------
                when S_OUTPUT =>
                    osum := ola_curr(to_integer(out_idx))
                          + ola_prev(to_integer(out_idx) + HOP_SIZE);

                    if osum > to_signed(2**(WIDTH-1)-1, WIDTH+1) then
                        pcm_out <= std_logic_vector(
                                     to_signed(2**(WIDTH-1)-1, WIDTH));
                    elsif osum < to_signed(-(2**(WIDTH-1)), WIDTH+1) then
                        pcm_out <= std_logic_vector(
                                     to_signed(-(2**(WIDTH-1)), WIDTH));
                    else
                        pcm_out <= std_logic_vector(osum(WIDTH-1 downto 0));
                    end if;
                    oval <= '1';

                    if out_idx = to_unsigned(HOP_SIZE - 1, 7) then
                        ola_prev <= ola_curr;
                        state    <= S_COLLECT;
                    else
                        out_idx <= out_idx + 1;
                    end if;

                when others =>
                    state <= S_COLLECT;

                end case;
            end if;
        end if;
    end process;

    pcm_out_valid <= oval;

end RTL;