--------------------------------------------------------------------------------
-- Precision_MAC_lut.vhd              (TCAS-II submission - cost model sweep)
--
-- Identical to Precision_MAC.vhd but with the use_dsp="no" synthesis
-- attribute applied to the product signal.  This forces Vivado to implement
-- the multiplier entirely in LUT fabric, preventing DSP48E1 inference.
--
-- PURPOSE:
--   Used exclusively for the J(B) hardware cost model characterisation sweep
--   (cost_model_sweep_v4.tcl).  The LUT count as a function of B provides a
--   technology-agnostic proxy for multiplier complexity that approximates
--   ASIC logic-element (LE) scaling.  Because array multipliers scale as
--   O(B²) in both LUT fabric and ASIC standard cells, the quadratic cost
--   model derived here generalises beyond the specific FPGA device.
--
--   The operational filterbank (MAC_Path_Wrapper, WOLA_Filterbank) uses the
--   standard Precision_MAC.vhd which ALLOWS DSP48 inference for maximum
--   efficiency.  These are two separate files serving two separate purposes:
--     Precision_MAC.vhd      → production use, DSP48 inference enabled
--     Precision_MAC_lut.vhd  → cost model sweep only, LUT-only forced
--
-- ATTRIBUTE:
--   use_dsp = "no" on the product signal prevents Vivado from absorbing
--   the multiply into a DSP48E1 primitive.  The synthesiser must implement
--   the full partial-product matrix in LUT6 slices, giving a clean LUT
--   count that scales quadratically with B as expected.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity Precision_MAC_lut is
    generic (
        WIDTH : integer := 24
    );
    port (
        clk     : in  std_logic;
        rst     : in  std_logic;
        en      : in  std_logic;
        a_in    : in  std_logic_vector(WIDTH-1 downto 0);
        b_in    : in  std_logic_vector(WIDTH-1 downto 0);
        res_out : out std_logic_vector(2*WIDTH-1 downto 0)
    );
end Precision_MAC_lut;

architecture RTL of Precision_MAC_lut is

    -- Prevent DSP48E1 inference so Vivado implements the full partial-product
    -- matrix in LUT fabric.  This is intentional - see file header comment.
    attribute use_dsp : string;
    signal product : signed(2*WIDTH-1 downto 0) := (others => '0');
    attribute use_dsp of product : signal is "no";

begin
    process(clk)
    begin
        if rising_edge(clk) then
            if rst = '1' then
                product <= (others => '0');
            elsif en = '1' then
                product <= signed(a_in) * signed(b_in);
            end if;
        end if;
    end process;
    res_out <= std_logic_vector(product);
end RTL;
