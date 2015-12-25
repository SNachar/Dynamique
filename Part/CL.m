function [M,K,F,U,Up,Upp]=CL(M,K,F,U,Up,Upp)
M(1,:)=[];
M(:,1)=[];
K(1,:)=[];
K(:,1)=[];
F(1,:)=[];
U(1,:)=[];
Up(1,:)=[];
Upp(1,:)=[];
Upp(:,1)=M\F(:,1);
end