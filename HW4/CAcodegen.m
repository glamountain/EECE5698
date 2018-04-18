function [ca_code]=CAcodegen(svnum)
% [ca_used]=CAcodegen(svnum), generates GPS L1 C/A Gold codes
%
% ca_used : a vector containing the desired output sequence
% svnum: Satellite number (can be a vector!), can generate from SV1 to
% SV32

k=length(svnum);            % number of satellites
ca_code=zeros(k,1023);
G1=zeros(1,1023);
G2=zeros(1,1023);

for i=1:k

  % the g2s vector holds the appropriate shift of the G2 code to
  % generate the C/A code (ex. for SV#19 -use a G2 shift of
  % g2s(19)=471) 
  
  g2s = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;469;470;471;...
	 472;473;474;509;512;513;514;515;516; 859;860;861;862]; 
  g2shift=g2s(svnum(i),1); 
  % ***** Generate G1 code ***** % load shift register 
  reg= -1*ones(1,10); 
  for n= 1:1023
    G1(n) = reg(10) ; 
    save1 = reg(3)*reg(10); 			
    reg(1,2:10) = reg(1:1:9) ;
    reg(1) = save1; 
  end
  % ***** Generate G2 code ***** % load shift register 
  reg= -1*ones(1,10);
  for n= 1:1023
    ????????
  end 
  % ***** Shift G2 code to get G2i ***** 
  G2i = ????????
  % ***** Form single sample C/Acode by multiplying G1 and G2 
  ca_code(i,:)=???????? 
end
