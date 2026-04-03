--------------------------------------------------------------------------------
-- MAC_Path_Wrapper.vhd
-- Synthesis wrapper for SAIF power analysis.
-- Instantiates only the precision-controlled MAC datapath.
-- Set this as the synthesis top in Vivado.
--------------------------------------------------------------------------------
library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity MAC_Path_Wrapper is
    port (
        clk         : in  std_logic;
        rst         : in  std_logic;
        -- Allocation config
        alloc_we    : in  std_logic;
        alloc_addr  : in  std_logic_vector(5 downto 0);
        alloc_data  : in  std_logic_vector(4 downto 0);
        -- Frame handshake
        frame_valid : in  std_logic;
        frame_done  : out std_logic;
        band_index  : out std_logic_vector(5 downto 0);
        band_valid  : out std_logic;
        band_ack    : in  std_logic;
        -- Operand inputs (from external stimulus)
        a_raw       : in  std_logic_vector(23 downto 0);
        b_raw       : in  std_logic_vector(23 downto 0);
        -- MAC output
        mac_result  : out std_logic_vector(47 downto 0)
    );
end MAC_Path_Wrapper;

architecture Structural of MAC_Path_Wrapper is
    signal mask      : std_logic_vector(23 downto 0);
    signal mac_en    : std_logic;
    signal a_masked  : std_logic_vector(23 downto 0);
    signal b_masked  : std_logic_vector(23 downto 0);
begin

    U_PC : entity work.Precision_Controller
        generic map (WIDTH => 24, NUM_BANDS => 64)
        port map (
            clk => clk, rst => rst,
            alloc_we => alloc_we, alloc_addr => alloc_addr,
            alloc_data => alloc_data,
            frame_valid => frame_valid, frame_done => frame_done,
            band_index => band_index, mask_out => mask,
            mac_en => mac_en, band_valid => band_valid,
            band_ack => band_ack);

    U_ISO : entity work.Operand_Isolation
        generic map (WIDTH => 24)
        port map (
            a_raw => a_raw, b_raw => b_raw,
            mask => mask,
            a_masked => a_masked, b_masked => b_masked);

    U_MAC : entity work.Precision_MAC
        generic map (WIDTH => 24)
        port map (
            clk => clk, rst => rst, en => mac_en,
            a_in => a_masked, b_in => b_masked,
            res_out => mac_result);

end Structural;