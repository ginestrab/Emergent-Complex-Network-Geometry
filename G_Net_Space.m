function [ a,R,k,t,chi] = G_Net_Space(N,m,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed by the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% If you use this code, please cite 
%  Z. Wu, G. Menichetti, C. Rahmede, G. Bianconi
% "Emergent Complex Network Geometry"
%   Scientific Reports 5, 10073 (2015).
% (c) Zhihao Wu, Giulia Menichetti, Christoph Rahmede and Ginestra Bianconi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code uses m=1,2,3,4...
% In the case in which m=1 the program returns the model for m=infinity.
%N maximal number of nodes
% a adjacency matrix of size N by N
%R vector of size N (curvature R(i) for each node  i
%k vector of size N vector of degrees k(i) of node i
%t vector of size N vector of number of  triangles t(i) for each node
% chi Euler characteristics when the network is of size N
					 
					 
	a=sparse(N,N);
    flag=1;
	t=zeros(1,N);
	a_occ=sparse(N,N);
	L=0;
	idaus=0;
	if m==1,
	 flag=0;
	 m=2;
	end
	for i1=1:3,
		for i2=i1+1:3,     
			L=L+1;
			a(i1,i2)=1;
			a(i2,i1)=1;
			a_occ(i1,i2)=m-1;  
			a_occ(i2,i1)=m-1;  
					 
		end
		t(i1)=1;
					 
	end
					 
	Lp=0;
	laus=4;
	for in=3+1:N,
		[I,J,V]=find(tril(a.*(a_occ>0)));
		norm=sum(V);
		x=rand(1)*norm;
		for nj1=1:max(size(V)),
			x=x-V(nj1);
			if x<0,
				nj=nj1;
			break;
		    end
		end
		l1=I(nj);
		l2=J(nj);
		t(l1)=t(l1)+1;
		t(l2)=t(l2)+1;
		t(in)=t(in)+1;
					 
		a_occ(l1,l2)=a_occ(l1,l2)-flag;
		a_occ(l2,l1)=a_occ(l2,l1)-flag;
					 
		L=L+1;
					 
		a(in,l1)=1;
		a(l1,in)=1;
		a_occ(in,l1)=m-1;  
		a_occ(l1,in)=m-1;  
					 
		L=L+1;
					 
		a(in,l2)=1;    
		a(l2,in)=1;
		a_occ(in,l2)=m-1;  
		a_occ(l2,in)=m-1;  
					 
		if ((rand(1)<p)),
			[I,J,V]=find(tril(a.*(a_occ>0)));
			norm=sum(V);
			x=rand(1)*norm;
			for nj1=1:max(size(V)),
				x=x-V(nj1);
				if x<0,
					 nj=nj1;
					 break;
				end
			end
			l1=I(nj);
			l2=J(nj);
			b=ones(N,1);
			b(l1,1)=0;
			b(l2,1)=0;
			b1=(a_occ(:,l1)>0);
			b2=(a_occ(:,l2)>0);
					 
			[J,V]=find(b.*(b1+b2-2*b1.*b2));
			if nnz(J)>0,
				ni=ceil(rand(1)*max(size(V)));
				i3=J(ni);
				if(a_occ(l1,i3)>0)
					 i1p=l2;
					 i2p=l1;
					 i3p=i3;
				end
				if(a_occ(l2,i3)>0)
					 i1p=l1;
					 i2p=l2;
					 i3p=i3;
				end
				B=a^2;
				Baus=(a(i1p,:).*(a_occ(i1p,:)>0))*(a(:,i3p).*(a_occ(:,i3p)>0));
				if ((a(i1p,i3p)==0) &&(Baus==B(i1p,i3p)))
					 
					 [J1,V1]=find(a(:,i1p).*a(:,i3p));
					 for n3=1:numel(V1),
						i2p=J1(n3);
						t(i1p)=t(i1p)+1;
						t(i2p)=t(i2p)+1;
						t(i3p)=t(i3p)+1;
						a_occ(i1p,i2p)=a_occ(i1p,i2p)-flag;
						a_occ(i2p,i1p)=a_occ(i2p,i1p)-flag;  
						a_occ(i3p,i2p)=a_occ(i3p,i2p)-flag;
						a_occ(i2p,i3p)=a_occ(i2p,i3p)-flag;  
					 end
					 
					 L=L+1; 
					 
					 a(i1p,i3p)=1;
					 a(i3p,i1p)=1;
					 a_occ(i1p,i3p)=m-1;  
					 a_occ(i3p,i1p)=m-1;
				end
			end
		end
					 
	end
					 
					 
	k=sum(a);
	R=(k>0)-k/2+t/3;
	chi=sum(R);

end
