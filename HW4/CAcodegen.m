function [ca_code] = CAcodegen(svnum)
% [ca_used]=CAcodegen(svnum), generates GPS L1 C/A Gold codes
%
% ca_used : a vector containing the desired output sequence
% svnum: Satellite number (can be a vector!), can generate from SV1 to
% SV32

k       = length(svnum);            % number of satellites
ca_code = zeros(k,1023);
G1      = zeros(1,1023);
G2      = zeros(1,1023);

for i=1:k
    
    % the g2s vector holds the appropriate shift of the G2 code to
    % generate the C/A code (ex. for SV#19 -use a G2 shift of
    % g2s(19)=471)
    
    g2s = [5;6;7;8;17;18;139;140;...
        141;251;252;254;255;256;257;258;...
        469;470;471;472;473;474;509;512;...
        513;514;515;516;859;860;861;862];
    g2shift = g2s(svnum(i),1);
    
    g1poly = [3;10];
    g2poly = [2;3;6;8;9;10];
    
    
    %% ***** Generate G1 and G2 codes *****
    g1_reg = -1 * ones(1,10);
    g2_reg = -1 * ones(1,10);
    for n= 1:1023
        
        % Load output register
        G1(n)       = g1_reg(10);
        G2(n)       = g2_reg(10);
        
        % Compute feedback and load shift register
        g1_reg(end) = prod(g1_reg(g1poly));
        g2_reg(end) = prod(g2_reg(g2poly));
        
        % Clock shift register
        g1_reg      = circshift(g1_reg,1,2);
        g2_reg      = circshift(g2_reg,1,2);
        
    end
    
    %% ***** Shift G2 code to get G2i *****
    G2i = circshift(G2,g2shift,2);
    
    
    %% ***** Form single sample C/Acode by multiplying G1 and G2
    ca_code(i,:) = G1 .* G2i;
    
end
