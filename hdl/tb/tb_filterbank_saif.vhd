--------------------------------------------------------------------------------
-- tb_filterbank_saif.vhd  (v2 - updated for mac_saif_sweep_v2.tcl)
--
-- G_RUN is now an entity generic so Vivado can override it via
-- -generic_top G_RUN=N without modifying the source file between runs.
--
-- Run index map (12 runs):
--   0  Uniform    40%   sum=410
--   1  Linear     40%   sum=410
--   2  Quadratic  40%   sum=410
--   3  Uniform    60%   sum=614
--   4  Linear     60%   sum=614
--   5  Quadratic  60%   sum=614
--   6  Uniform    80%   sum=819
--   7  Linear     80%   sum=819
--   8  Quadratic  80%   sum=819
--   9  Uniform   100%   sum=1024
--  10  Linear    100%   sum=1024
--  11  Quadratic 100%   sum=1024
--
-- Allocation vectors sourced from E1_iso_bit_allocation_<N>pct.csv
-- (ALPHA=1.37, BETA=0 cost model, updated run).
-- All sums verified: 40%->410, 60%->614, 80%->819, 100%->1024.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity tb_filterbank_saif is
    generic (
        G_RUN : integer := 0   -- overridden by Vivado -generic_top G_RUN=N
    );
end entity tb_filterbank_saif;

architecture sim of tb_filterbank_saif is

    constant CLK_PERIOD : time    := 10 ns;
    constant NUM_BANDS  : integer := 64;
    constant NUM_FRAMES : integer := 256;

    -- 12 runs x 64 bands
    type alloc_table_t is array (0 to 11, 0 to 63) of integer range 1 to 24;
    constant ALLOC : alloc_table_t := (

        -- Run 0: Uniform 40%  sum=410
        (7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
         6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
         6,6,6,6,6,6,6,6,6,6,6,6),

        -- Run 1: Linear 40%  sum=410
        (1,4,6,7,8,9,8,8,8,8,8,8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,
         7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,
         5,5,5,5,5,5,5,5,5,5,5,5),

        -- Run 2: Quadratic 40%  sum=410
        (1,4,6,7,9,9,9,9,9,9,8,8,9,8,8,8,7,7,7,7,7,7,7,7,7,7,
         7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,5,5,5,
         5,5,5,5,5,5,5,5,5,5,5,5),

        -- Run 3: Uniform 60%  sum=614
        (10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
         9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9),

        -- Run 4: Linear 60%  sum=614
        (2,7,9,10,11,13,13,12,12,11,11,11,11,11,11,11,11,11,11,11,
         11,11,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
         9,9,9,9,9,9,9,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8),

        -- Run 5: Quadratic 60%  sum=614
        (3,7,9,10,12,12,12,12,12,12,11,12,12,11,11,11,11,11,11,11,
         11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,9,9,9,
         9,9,9,9,9,9,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8),

        -- Run 6: Uniform 80%  sum=819
        (13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
         13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
         13,13,13,13,13,13,13,13,13,13,13,
         12,12,12,12,12,12,12,12,12,12,12,12,12),

        -- Run 7: Linear 80%  sum=819
        (5,10,13,14,15,16,16,16,15,14,14,14,14,14,14,14,14,14,14,14,
         14,14,14,14,14,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
         12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
         11,11,11,11,11,11,11,11),

        -- Run 8: Quadratic 80%  sum=819
        (6,10,13,14,15,16,16,15,15,15,14,15,15,14,14,14,14,14,14,14,
         14,14,14,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
         12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
         11,11,11,11,11,11,11,11,11),

        -- Run 9: Uniform 100%  sum=1024
        (16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
         16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
         16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
         16,16,16,16),

        -- Run 10: Linear 100%  sum=1024
        (8,14,16,17,18,17,17,17,18,18,18,18,18,18,18,17,17,17,17,17,
         17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,16,16,
         16,16,16,16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
         14,14,14,14,14),

        -- Run 11: Quadratic 100%  sum=1024
        (9,14,16,17,18,19,18,18,18,18,17,17,18,17,17,17,17,17,17,17,
         17,17,17,17,17,17,17,17,17,17,16,16,16,16,16,16,16,16,16,16,
         16,16,16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
         14,14,14,14,14)
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

    type label_array_t is array (0 to 11) of string(1 to 14);
    constant RUN_LABEL : label_array_t := (
        "Uniform    40%",
        "Linear     40%",
        "Quadratic  40%",
        "Uniform    60%",
        "Linear     60%",
        "Quadratic  60%",
        "Uniform    80%",
        "Linear     80%",
        "Quadratic  80%",
        "Uniform   100%",
        "Linear    100%",
        "Quadratic 100%"
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