--------------------------------------------------------------------------------
-- Precision_Controller.vhd              (TCAS-II submission - full filterbank)
--
-- Four-state FSM + 64×5-bit distributed-RAM LUT.
--
-- Per Section IV of the paper, the FSM sequences through all 64 subbands
-- once per frame.  For each band k it:
--   S_IDLE    → waits for frame_valid
--   S_FETCH   → reads B_k from the allocation LUT (1-cycle read latency)
--   S_EXECUTE → drives mask_out and mac_en for one clock
--   S_WAIT    → holds until band_ack from the gain sub-FSM, then increments k
--
-- Mask generation:  mask(i) = '1' for i ∈ [WIDTH-1 : WIDTH-B_k],
--                              '0' otherwise.
-- Implemented as a single barrel-shift of all-ones.
-- Area: ~48 LUTs on Artix-7 (matches paper Table II).
--
-- The allocation vector is loaded at startup through a simple synchronous
-- write port and remains static during normal operation.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity Precision_Controller is
    generic (
        WIDTH     : integer := 24;
        NUM_BANDS : integer := 64
    );
    port (
        clk         : in  std_logic;
        rst         : in  std_logic;
        -- Configuration write port
        alloc_we    : in  std_logic;
        alloc_addr  : in  std_logic_vector(5 downto 0);
        alloc_data  : in  std_logic_vector(4 downto 0);
        -- Frame handshake
        frame_valid : in  std_logic;
        frame_done  : out std_logic;
        -- Per-band control (active during S_EXECUTE)
        band_index  : out std_logic_vector(5 downto 0);
        mask_out    : out std_logic_vector(WIDTH-1 downto 0);
        mac_en      : out std_logic;
        band_valid  : out std_logic;
        -- Back-pressure from gain sub-FSM
        band_ack    : in  std_logic
    );
end Precision_Controller;

architecture RTL of Precision_Controller is

    type state_t is (S_IDLE, S_FETCH, S_EXECUTE, S_WAIT);
    signal state : state_t := S_IDLE;

    signal k_cnt : unsigned(5 downto 0) := (others => '0');

    type lut_t is array (0 to NUM_BANDS-1) of unsigned(4 downto 0);
    signal alloc_lut : lut_t := (others => to_unsigned(16, 5));

    signal bk_reg : unsigned(4 downto 0) := (others => '0');
    signal mask_i : std_logic_vector(WIDTH-1 downto 0);

begin

    -- Allocation LUT write -----------------------------------------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if alloc_we = '1' then
                alloc_lut(to_integer(unsigned(alloc_addr))) <=
                    unsigned(alloc_data);
            end if;
        end if;
    end process;

    -- Mask generation (combinational) ------------------------------------------
    process(bk_reg)
        variable n : integer range 0 to WIDTH;
    begin
        n := to_integer(bk_reg);
        if n = 0 then
            mask_i <= (others => '0');
        elsif n >= WIDTH then
            mask_i <= (others => '1');
        else
            mask_i <= std_logic_vector(
                shift_left(to_unsigned(2**WIDTH - 1, WIDTH), WIDTH - n));
        end if;
    end process;

    -- Four-state FSM (S_IDLE, S_FETCH, S_EXECUTE, S_WAIT) ----------------------
    process(clk)
    begin
        if rising_edge(clk) then
            if rst = '1' then
                state      <= S_IDLE;
                k_cnt      <= (others => '0');
                bk_reg     <= (others => '0');
                mac_en     <= '0';
                band_valid <= '0';
                frame_done <= '0';
            else
                mac_en     <= '0';
                band_valid <= '0';
                frame_done <= '0';

                case state is
                when S_IDLE =>
                    k_cnt <= (others => '0');
                    if frame_valid = '1' then
                        state <= S_FETCH;
                    end if;

                when S_FETCH =>
                    bk_reg <= alloc_lut(to_integer(k_cnt));
                    state  <= S_EXECUTE;

                when S_EXECUTE =>
                    mac_en     <= '1';
                    band_valid <= '1';
                    state      <= S_WAIT;

                when S_WAIT =>
                    -- Hold until the gain sub-FSM acknowledges completion
                    if band_ack = '1' then
                        if k_cnt = to_unsigned(NUM_BANDS-1, 6) then
                            frame_done <= '1';
                            state      <= S_IDLE;
                        else
                            k_cnt <= k_cnt + 1;
                            state <= S_FETCH;
                        end if;
                    end if;
                end case;
            end if;
        end if;
    end process;

    band_index <= std_logic_vector(k_cnt);
    mask_out   <= mask_i;

end RTL;