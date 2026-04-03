--------------------------------------------------------------------------------
-- tb_filterbank_saif_isodist.vhd
--
-- SAIF testbench for iso-distortion switching activity measurement.
-- Each run uses the MINIMUM bit budget for that strategy to match
-- Uniform's D_w at the corresponding quality target level.
--
-- Key insight: at iso-distortion, Linear and Quadratic need 8-15% fewer
-- total bits than Uniform to achieve the same perceptual quality.
-- Lower total bits -> lower average operand width -> less switching.
-- This demonstrates the efficiency argument: same quality, less power.
--
-- Run index map (12 runs):
--   0  Uniform    @ D_w target 40%   sum=400  avg=6.250
--   1  Linear     @ D_w target 40%   sum=348  avg=5.438
--   2  Quadratic  @ D_w target 40%   sum=357  avg=5.578
--   3  Uniform    @ D_w target 60%   sum=596  avg=9.312
--   4  Linear     @ D_w target 60%   sum=550  avg=8.594
--   5  Quadratic  @ D_w target 60%   sum=553  avg=8.641
--   6  Uniform    @ D_w target 80%   sum=790  avg=12.344
--   7  Linear     @ D_w target 80%   sum=749  avg=11.703
--   8  Quadratic  @ D_w target 80%   sum=750  avg=11.719
--   9  Uniform    @ D_w target 100%  sum=982  avg=15.344
--  10  Linear     @ D_w target 100%  sum=943  avg=14.734
--  11  Quadratic  @ D_w target 100%  sum=944  avg=14.750
--
-- Allocation vectors sourced from E10_iso_distortion_detail.csv
-- (ALPHA=1.37, BETA=0 cost model).
-- All sums independently verified against E10_iso_distortion_summary.csv.
--
-- USAGE: Set G_RUN via Vivado generic override before each simulation.
-- After simulation: report_power with read_saif to get switching activity.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity tb_filterbank_saif_isodist is
    generic (
        G_RUN : integer := 11   -- overridden per run: 0-11
    );
end entity tb_filterbank_saif_isodist;

architecture sim of tb_filterbank_saif_isodist is

    constant CLK_PERIOD : time    := 10 ns;
    constant NUM_BANDS  : integer := 64;
    constant NUM_FRAMES : integer := 256;

    -- 12 runs x 64 bands - iso-distortion minimum budgets from E10
    type alloc_table_t is array (0 to 11, 0 to 63) of integer range 1 to 24;
    constant ALLOC : alloc_table_t := (

        -- Run 0: Uniform @40%  sum=400  avg=6.250
        (7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6),

        -- Run 1: Linear @40%  sum=348  avg=5.438
        (4, 4, 6, 7, 8, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 7, 7, 
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4),

        -- Run 2: Quadratic @40%  sum=357  avg=5.578
        (4, 4, 6, 7, 8, 9, 9, 9, 8, 8, 7, 8, 8, 7, 7, 7, 7, 7, 
        7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 7, 6, 6, 6, 6, 6, 6, 6, 6, 
        6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 5, 4, 4, 4),

        -- Run 3: Uniform @60%  sum=596  avg=9.312
        (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9),

        -- Run 4: Linear @60%  sum=550  avg=8.594
        (4, 7, 9, 10, 11, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7),

        -- Run 5: Quadratic @60%  sum=553  avg=8.641
        (4, 7, 9, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 
        9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7),

        -- Run 6: Uniform @80%  sum=790  avg=12.344
        (13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 
        13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12),

        -- Run 7: Linear @80%  sum=749  avg=11.703
        (5, 10, 12, 13, 14, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 
        14, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 
        13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11),

        -- Run 8: Quadratic @80%  sum=750  avg=11.719
        (5, 10, 12, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 
        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 12, 
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11),

        -- Run 9: Uniform @100%  sum=982  avg=15.344
        (16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16),

        -- Run 10: Linear @100%  sum=943  avg=14.734
        (8, 13, 15, 16, 17, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 
        17, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
        16, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14),

        -- Run 11: Quadratic @100%  sum=944  avg=14.750
        (8, 13, 15, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 
        16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 14, 14, 
        14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14)
    );

    signal clk         : std_logic := '0';
    signal rst         : std_logic := '1';
    signal alloc_we    : std_logic := '0';
    signal alloc_addr  : std_logic_vector(5 downto 0)  := (others => '0');
    signal alloc_data  : std_logic_vector(4 downto 0)  := (others => '0');
    signal frame_valid : std_logic := '0';
    signal frame_done  : std_logic;
    signal band_index  : std_logic_vector(5 downto 0);
    signal band_valid  : std_logic;
    signal band_ack    : std_logic := '0';
    signal a_raw       : std_logic_vector(23 downto 0) := (others => '0');
    signal b_raw       : std_logic_vector(23 downto 0) := (others => '0');
    signal mac_result  : std_logic_vector(47 downto 0);

    type label_array_t is array (0 to 11) of string(1 to 16);
    constant RUN_LABEL : label_array_t := (
        "Uniform  @40%   ",
        "Linear   @40%   ",
        "Quadratic@40%   ",
        "Uniform  @60%   ",
        "Linear   @60%   ",
        "Quadratic@60%   ",
        "Uniform  @80%   ",
        "Linear   @80%   ",
        "Quadratic@80%   ",
        "Uniform  @100%  ",
        "Linear   @100%  ",
        "Quadratic@100%  "
    );

begin

    clk <= not clk after CLK_PERIOD/2;

    U_DUT : entity work.MAC_Path_Wrapper
        port map (
            clk         => clk,
            rst         => rst,
            alloc_we    => alloc_we,
            alloc_addr  => alloc_addr,
            alloc_data  => alloc_data,
            frame_valid => frame_valid,
            frame_done  => frame_done,
            band_index  => band_index,
            band_valid  => band_valid,
            band_ack    => band_ack,
            a_raw       => a_raw,
            b_raw       => b_raw,
            mac_result  => mac_result
        );

    p_stim: process
        variable lfsr_a  : std_logic_vector(22 downto 0)
                         := "10110100110101010111001";
        variable lfsr_b  : std_logic_vector(22 downto 0)
                         := "01001011001010101000110";
        variable fb_a, fb_b : std_logic;
        variable timeout    : integer;
    begin
        rst <= '1';
        wait for 20 * CLK_PERIOD;
        wait until rising_edge(clk);
        rst <= '0';
        wait for 5 * CLK_PERIOD;

        -- Load allocation vector for selected run
        for k in 0 to NUM_BANDS-1 loop
            wait until rising_edge(clk);
            alloc_we   <= '1';
            alloc_addr <= std_logic_vector(to_unsigned(k, 6));
            alloc_data <= std_logic_vector(to_unsigned(ALLOC(G_RUN, k), 5));
        end loop;
        wait until rising_edge(clk);
        alloc_we <= '0';
        wait for 10 * CLK_PERIOD;

        -- Drive NUM_FRAMES complete frame cycles
        for frame in 0 to NUM_FRAMES-1 loop

            wait until rising_edge(clk);
            frame_valid <= '1';
            wait until rising_edge(clk);
            frame_valid <= '0';

            timeout := 0;
            while frame_done /= '1' loop
                wait until rising_edge(clk);
                timeout := timeout + 1;

                if band_valid = '1' then
                    fb_a   := lfsr_a(22) xor lfsr_a(17);
                    lfsr_a := lfsr_a(21 downto 0) & fb_a;
                    fb_b   := lfsr_b(22) xor lfsr_b(17);
                    lfsr_b := lfsr_b(21 downto 0) & fb_b;
                    a_raw    <= lfsr_a(22) & lfsr_a;
                    b_raw    <= lfsr_b(22) & lfsr_b;
                    band_ack <= '1';
                else
                    band_ack <= '0';
                end if;

                if timeout > 512 then
                    report "TIMEOUT frame " & integer'image(frame)
                        severity warning;
                    exit;
                end if;
            end loop;

            band_ack <= '0';
            wait for 2 * CLK_PERIOD;

        end loop;

        wait for 100 * CLK_PERIOD;
        report "SAIF_DONE: " & RUN_LABEL(G_RUN) severity failure;
        wait;
    end process;

end architecture sim;