--------------------------------------------------------------------------------
-- Operand_Isolation.vhd                  (TCAS-II submission - full filterbank)
--
-- Combinational AND-gate arrays for per-band operand isolation.
--
-- For each MAC operand, every bit is ANDed with the corresponding mask bit
-- produced by the Precision Controller.  Bits below the allocated precision
-- B_k are forced to logic-0, eliminating all switching activity in the
-- corresponding partial-product rows of the downstream DSP48 multiplier.
--
-- Two independent WIDTH-bit AND arrays (one per operand).
-- Purely combinational - no registers, no clock.
--
-- Area: 2 × WIDTH LUT1s ≈ 24 LUTs in the packed Artix-7 fabric.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity Operand_Isolation is
    generic (
        WIDTH : integer := 24
    );
    port (
        a_raw    : in  std_logic_vector(WIDTH-1 downto 0);
        b_raw    : in  std_logic_vector(WIDTH-1 downto 0);
        mask     : in  std_logic_vector(WIDTH-1 downto 0);
        a_masked : out std_logic_vector(WIDTH-1 downto 0);
        b_masked : out std_logic_vector(WIDTH-1 downto 0)
    );
end Operand_Isolation;

architecture Combinational of Operand_Isolation is
begin
    a_masked <= a_raw and mask;
    b_masked <= b_raw and mask;
end Combinational;