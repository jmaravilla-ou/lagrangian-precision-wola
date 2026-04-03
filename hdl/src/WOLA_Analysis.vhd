--------------------------------------------------------------------------------
-- WOLA_Analysis_v2.vhd                    (TCAS-II submission - full filterbank)
--
-- Synthesisable WOLA analysis front-end using Xilinx xfft_0 IP (v9.1.15).
-- Replaces the behavioural real-arithmetic FFT in WOLA_Analysis.vhd.
--
-- PIPELINE:
--   Input buffer (circular) → Hann window → xfft_0 (128-pt FFT) → output
--
-- XFFT_0 INTERFACE (from xfft_0.vhd):
--   s_axis_data_tdata  [47:0]  : {Re[23:0], Im[23:0]} packed
--                                Im in [23:0], Re in [47:24]
--   s_axis_data_tvalid         : assert to push one sample
--   s_axis_data_tready         : IP ready to accept (backpressure)
--   s_axis_data_tlast          : assert on sample 127 (last of frame)
--   m_axis_data_tdata  [47:0]  : output, same packing
--   m_axis_data_tvalid         : output sample valid
--   m_axis_data_tready         : we keep this '1' (always accept)
--   m_axis_data_tlast          : pulses on last output sample (bin 127)
--   s_axis_config_tdata[15:0]  : bit[0]=FWD_INV (0=FFT), bits[8:1]=SCALE_SCH
--   s_axis_config_tvalid       : send config once at startup
--   s_axis_config_tready       : IP ready to accept config
--
-- SCALING:
--   128-point FFT has 7 stages. SCALE_SCH = 8'b00101010 (bits[8:1] = 0x2A)
--   giving /2 at stages 2,4,6 → total right-shift of 3 (divide by 8).
--   This prevents overflow on speech-level inputs with 24-bit data.
--   The synthesis IFFT uses complementary scaling (none) to restore level.
--
-- NATURAL ORDER:
--   C_HAS_NATURAL_INPUT=1, C_HAS_NATURAL_OUTPUT=1 → samples fed in order
--   0,1,...,127; bins received in order 0,1,...,127. No bit-reversal needed.
--
-- OUTPUT:
--   We emit bins 0..63 only (one-sided spectrum of real signal).
--   sub_valid pulses once per bin. frame_complete pulses after bin 63.
--   sub_mag is an alpha-max-beta-min approximation of magnitude.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use IEEE.MATH_REAL.ALL;

entity WOLA_Analysis is
    generic (
        WIDTH      : integer := 24;
        FRAME_LEN  : integer := 128;
        HOP_SIZE   : integer := 64;
        NUM_BANDS  : integer := 64
    );
    port (
        clk            : in  std_logic;
        rst            : in  std_logic;
        pcm_in         : in  std_logic_vector(WIDTH-1 downto 0);
        pcm_in_valid   : in  std_logic;
        sub_re         : out std_logic_vector(WIDTH-1 downto 0);
        sub_im         : out std_logic_vector(WIDTH-1 downto 0);
        sub_mag        : out std_logic_vector(WIDTH-1 downto 0);
        sub_index      : out std_logic_vector(5 downto 0);
        sub_valid      : out std_logic;
        frame_valid    : out std_logic;     -- pulses on bin 0 of each frame
        frame_complete : out std_logic      -- pulses after bin 63
    );
end WOLA_Analysis;

architecture RTL of WOLA_Analysis is

    ---------------------------------------------------------------------------
    -- xfft_0 component declaration (matches xfft_0.vhd entity exactly)
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
    -- Hann window ROM
    -- w[n] = 0.5*(1 - cos(2*pi*n/128)), quantised to WIDTH-1 bits (unsigned)
    ---------------------------------------------------------------------------
    type rom_t is array (0 to FRAME_LEN-1) of signed(WIDTH-1 downto 0);
    function init_hann return rom_t is
        variable r : rom_t;
    begin
        for n in 0 to FRAME_LEN-1 loop
            r(n) := to_signed(integer(round(
                0.5*(1.0 - cos(2.0*MATH_PI*real(n)/real(FRAME_LEN)))
                * real(2**(WIDTH-1) - 1))), WIDTH);
        end loop;
        return r;
    end function;
    constant HANN_ROM : rom_t := init_hann;

    ---------------------------------------------------------------------------
    -- Circular sample buffer (FRAME_LEN deep, WIDTH wide)
    ---------------------------------------------------------------------------
    type sbuf_t is array (0 to FRAME_LEN-1) of signed(WIDTH-1 downto 0);
    signal circ_buf     : sbuf_t := (others => (others => '0'));
    signal wr_ptr       : unsigned(6 downto 0) := (others => '0');
    signal hop_cnt      : unsigned(6 downto 0) := (others => '0');

    -- Handshake between input process and processing FSM
    signal frame_pending  : std_logic := '0';
    signal frame_ack      : std_logic := '0';
    signal snap_rd_base   : unsigned(6 downto 0) := (others => '0');

    ---------------------------------------------------------------------------
    -- Processing FSM states
    -- IDLE     → wait for frame_pending
    -- WINDOW   → read 128 samples, multiply by Hann, store to win_buf
    -- FEED_CFG → send config word to xfft_0 (FWD, scale schedule)
    -- FEED     → stream 128 windowed samples to xfft_0
    -- COLLECT  → accept 128 output bins from xfft_0, emit bins 0..63
    ---------------------------------------------------------------------------
    type state_t is (S_IDLE, S_WINDOW, S_FEED_CFG, S_FEED, S_COLLECT);
    signal state : state_t := S_IDLE;

    -- Windowed frame buffer
    type wbuf_t is array (0 to FRAME_LEN-1) of signed(WIDTH-1 downto 0);
    signal win_buf   : wbuf_t := (others => (others => '0'));
    signal win_idx   : unsigned(6 downto 0) := (others => '0');
    signal rd_base   : unsigned(6 downto 0) := (others => '0');

    -- Feed-to-xfft counters and signals
    signal feed_idx  : unsigned(6 downto 0) := (others => '0');
    signal cfg_sent  : std_logic := '0';      -- config sent flag (one-time)

    -- xfft_0 input signals
    signal fft_din_tdata  : std_logic_vector(47 downto 0) := (others => '0');
    signal fft_din_tvalid : std_logic := '0';
    signal fft_din_tready : std_logic;
    signal fft_din_tlast  : std_logic := '0';

    -- xfft_0 config signals
    signal fft_cfg_tdata  : std_logic_vector(15 downto 0) := (others => '0');
    signal fft_cfg_tvalid : std_logic := '0';
    signal fft_cfg_tready : std_logic;

    -- xfft_0 output signals
    signal fft_dout_tdata  : std_logic_vector(47 downto 0);
    signal fft_dout_tvalid : std_logic;
    signal fft_dout_tlast  : std_logic;

    -- Collect-from-xfft counter
    signal col_idx   : unsigned(6 downto 0) := (others => '0');

    -- Output registers
    signal o_re       : std_logic_vector(WIDTH-1 downto 0) := (others => '0');
    signal o_im       : std_logic_vector(WIDTH-1 downto 0) := (others => '0');
    signal o_mag      : std_logic_vector(WIDTH-1 downto 0) := (others => '0');
    signal o_idx      : std_logic_vector(5 downto 0)       := (others => '0');
    signal o_valid    : std_logic := '0';
    signal o_fvalid   : std_logic := '0';
    signal o_fcomplete: std_logic := '0';

begin

    ---------------------------------------------------------------------------
    -- xfft_0 instantiation
    -- Config word: bits[8:1] = SCALE_SCH = 8'b00101010 = 0x2A (scale /8)
    --              bit[0]    = FWD_INV   = 0 (forward FFT)
    -- Full 16-bit config = 0x0054  (0x2A << 1 | 0)
    ---------------------------------------------------------------------------
    fft_cfg_tdata <= x"0054";   -- FWD FFT, scale schedule 0x2A

    U_FFT : xfft_0
        port map (
            aclk                        => clk,
            s_axis_config_tdata         => fft_cfg_tdata,
            s_axis_config_tvalid        => fft_cfg_tvalid,
            s_axis_config_tready        => fft_cfg_tready,
            s_axis_data_tdata           => fft_din_tdata,
            s_axis_data_tvalid          => fft_din_tvalid,
            s_axis_data_tready          => fft_din_tready,
            s_axis_data_tlast           => fft_din_tlast,
            m_axis_data_tdata           => fft_dout_tdata,
            m_axis_data_tvalid          => fft_dout_tvalid,
            m_axis_data_tready          => '1',   -- always accept output
            m_axis_data_tlast           => fft_dout_tlast,
            event_frame_started         => open,
            event_tlast_unexpected      => open,
            event_tlast_missing         => open,
            event_status_channel_halt   => open,
            event_data_in_channel_halt  => open,
            event_data_out_channel_halt => open
        );

    ---------------------------------------------------------------------------
    -- Input buffer process
    -- Runs independently of processing FSM.
    -- Accepts one sample per clock when pcm_in_valid = '1'.
    -- Sets frame_pending every HOP_SIZE samples.
    ---------------------------------------------------------------------------
    p_input: process(clk)
    begin
        if rising_edge(clk) then
            if rst = '1' then
                wr_ptr        <= (others => '0');
                hop_cnt       <= (others => '0');
                frame_pending <= '0';
                snap_rd_base  <= (others => '0');
            else
                if frame_ack = '1' then
                    frame_pending <= '0';
                end if;

                if pcm_in_valid = '1' then
                    circ_buf(to_integer(wr_ptr)) <= signed(pcm_in);
                    wr_ptr <= wr_ptr + 1;

                    if hop_cnt = to_unsigned(HOP_SIZE - 1, 7) then
                        hop_cnt       <= (others => '0');
                        frame_pending <= '1';
                        snap_rd_base  <= wr_ptr + 1;
                    else
                        hop_cnt <= hop_cnt + 1;
                    end if;
                end if;
            end if;
        end if;
    end process;

    ---------------------------------------------------------------------------
    -- Processing FSM
    ---------------------------------------------------------------------------
    p_proc: process(clk)
        variable samp    : signed(WIDTH-1 downto 0);
        variable win_c   : signed(WIDTH-1 downto 0);
        variable prod    : signed(2*WIDTH-1 downto 0);
        variable rd_addr : unsigned(6 downto 0);
        -- Magnitude variables
        variable re_v, im_v         : signed(WIDTH-1 downto 0);
        variable abs_re_v, abs_im_v : signed(WIDTH-1 downto 0);
        variable abs_r, abs_i       : unsigned(WIDTH-2 downto 0);
        variable mx, mn             : unsigned(WIDTH-2 downto 0);
    begin
        if rising_edge(clk) then
            if rst = '1' then
                state         <= S_IDLE;
                win_idx       <= (others => '0');
                feed_idx      <= (others => '0');
                col_idx       <= (others => '0');
                cfg_sent      <= '0';
                frame_ack     <= '0';
                fft_din_tvalid<= '0';
                fft_din_tlast <= '0';
                fft_cfg_tvalid<= '0';
                o_valid       <= '0';
                o_fvalid      <= '0';
                o_fcomplete   <= '0';
            else
                -- Default pulse signals
                frame_ack   <= '0';
                o_valid     <= '0';
                o_fvalid    <= '0';
                o_fcomplete <= '0';

                case state is

                -- ----------------------------------------------------------------
                when S_IDLE =>
                    fft_din_tvalid <= '0';
                    fft_din_tlast  <= '0';
                    if frame_pending = '1' then
                        frame_ack <= '1';
                        rd_base   <= snap_rd_base;
                        win_idx   <= (others => '0');
                        state     <= S_WINDOW;
                    end if;

                -- ----------------------------------------------------------------
                -- S_WINDOW: one sample per clock, 128 clocks total
                -- Multiply sample by Hann coefficient, store to win_buf
                -- Product is Q1.23 × Q0.23 → take [2*WIDTH-2:WIDTH-1] = Q1.22
                -- ----------------------------------------------------------------
                when S_WINDOW =>
                    rd_addr := rd_base + win_idx;
                    samp    := circ_buf(to_integer(rd_addr));
                    win_c   := HANN_ROM(to_integer(win_idx));
                    prod    := samp * win_c;
                    win_buf(to_integer(win_idx)) <= prod(2*WIDTH-2 downto WIDTH-1);

                    if win_idx = to_unsigned(FRAME_LEN - 1, 7) then
                        -- Send config to xfft before feeding data
                        -- (only needed once; re-send each frame for safety)
                        fft_cfg_tvalid <= '1';
                        state <= S_FEED_CFG;
                    else
                        win_idx <= win_idx + 1;
                    end if;

                -- ----------------------------------------------------------------
                -- S_FEED_CFG: hold config valid until IP asserts ready
                -- ----------------------------------------------------------------
                when S_FEED_CFG =>
                    if fft_cfg_tready = '1' then
                        fft_cfg_tvalid <= '0';
                        feed_idx       <= (others => '0');
                        state          <= S_FEED;
                    end if;

                -- ----------------------------------------------------------------
                -- S_FEED: stream 128 windowed (real) samples to xfft_0
                -- Im channel = 0 (real-valued input)
                -- Data packing: [47:24]=Re, [23:0]=Im
                -- Advance only when tready is asserted (backpressure respected)
                -- ----------------------------------------------------------------
                when S_FEED =>
                    if fft_din_tready = '1' or fft_din_tvalid = '0' then
                        -- Pack: Im=0, Re=windowed sample
                        fft_din_tdata  <= std_logic_vector(
                            win_buf(to_integer(feed_idx))) & x"000000";
                        fft_din_tvalid <= '1';
                        fft_din_tlast  <= '0';

                        if feed_idx = to_unsigned(FRAME_LEN - 2, 7) then
                            -- Next cycle will be last sample
                            null;
                        end if;

                        if feed_idx = to_unsigned(FRAME_LEN - 1, 7) then
                            fft_din_tlast <= '1';
                            fft_din_tvalid<= '1';
                            col_idx       <= (others => '0');
                            state         <= S_COLLECT;
                        else
                            feed_idx <= feed_idx + 1;
                        end if;
                    end if;

                -- ----------------------------------------------------------------
                -- S_COLLECT: accept xfft_0 output bins 0..127
                -- Emit bins 0..63 to downstream (one-sided spectrum)
                -- Output packing from xfft: [47:24]=Re, [23:0]=Im
                -- ----------------------------------------------------------------
                when S_COLLECT =>
                    -- De-assert feed signals after last sample sent
                    fft_din_tvalid <= '0';
                    fft_din_tlast  <= '0';

                    if fft_dout_tvalid = '1' then
                        if col_idx < to_unsigned(NUM_BANDS, 7) then
                            -- Bins 0..63: emit to downstream
                            re_v := signed(fft_dout_tdata(47 downto 24));
                            im_v := signed(fft_dout_tdata(23 downto 0));

                            -- Alpha-max beta-min magnitude approximation
                            -- |z| ≈ max(|Re|,|Im|) + 3/8 * min(|Re|,|Im|)
                            abs_re_v := abs(re_v);
                            abs_im_v := abs(im_v);
                            abs_r    := unsigned(abs_re_v(WIDTH-2 downto 0));
                            abs_i    := unsigned(abs_im_v(WIDTH-2 downto 0));
                            if abs_r >= abs_i then
                                mx := abs_r; mn := abs_i;
                            else
                                mx := abs_i; mn := abs_r;
                            end if;

                            o_re    <= std_logic_vector(re_v);
                            o_im    <= std_logic_vector(im_v);
                            o_mag   <= std_logic_vector(
                                         signed('0' & std_logic_vector(
                                         mx + shift_right(mn,2)
                                            + shift_right(mn,3))));
                            o_idx   <= std_logic_vector(col_idx(5 downto 0));
                            o_valid <= '1';

                            if col_idx = 0 then
                                o_fvalid <= '1';   -- frame_valid on first bin
                            end if;

                            if col_idx = to_unsigned(NUM_BANDS - 1, 7) then
                                o_fcomplete <= '1';
                            end if;
                        end if;
                        -- Bins 64..127 are discarded (conjugate mirror)

                        if col_idx = to_unsigned(FRAME_LEN - 1, 7) then
                            state <= S_IDLE;
                        else
                            col_idx <= col_idx + 1;
                        end if;
                    end if;

                when others =>
                    state <= S_IDLE;

                end case;
            end if;
        end if;
    end process;

    -- Output assignments
    sub_re         <= o_re;
    sub_im         <= o_im;
    sub_mag        <= o_mag;
    sub_index      <= o_idx;
    sub_valid      <= o_valid;
    frame_valid    <= o_fvalid;
    frame_complete <= o_fcomplete;

end RTL;