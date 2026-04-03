--------------------------------------------------------------------------------
-- Precision_MAC.vhd                      (TCAS-II submission - full filterbank)
--
-- Signed multiply-accumulate unit (single-cycle multiply, registered output).
-- DSP48E1 inference is allowed (no use_dsp="no" attribute) so Vivado will
-- map the 24×24 signed multiply into one DSP48 slice.
--
-- The enable input (en) is driven by the Precision Controller's mac_en signal;
-- when de-asserted the output register holds its previous value and the
-- multiplier inputs should already be zero-masked by the Operand_Isolation
-- block, further suppressing switching activity.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity Precision_MAC is
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
end Precision_MAC;

architecture RTL of Precision_MAC is
    signal product : signed(2*WIDTH-1 downto 0) := (others => '0');
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