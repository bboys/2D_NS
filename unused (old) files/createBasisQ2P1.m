function [ localMatrix ] = createBasisQ2P1(  )
% creates the necessary matrices (mass, stiffness and 3rd order) for the 2D
% NS problem on a 2D mesh discretised with quadrilaterals both for the
% velocity unkowns (Q2) and pressure (P1)

nrVBasisF = 9;
nrPBasisF = 3;

% velocity

% 1D basis functions
syms t x y
points = [0,1/2,1]; % define the reference element ([0,1]^2)
lambda1 = (t-points(2))*(t-points(3)); 
lambda1 = lambda1/subs(lambda1,t,points(1));
lambda2 = (t-points(1))*(t-points(3)); 
lambda2 = lambda2/subs(lambda2,t,points(2));
lambda3 = (t-points(1))*(t-points(2)); 
lambda3 = lambda3/subs(lambda3,t,points(3));

%2D basis functions (labeled as:  bottom > top, left > right) (Q2)
basisF(1) = subs(lambda1,t,y)*subs(lambda1,t,x);
basisF(2) = subs(lambda2,t,y)*subs(lambda1,t,x);
basisF(3) = subs(lambda3,t,y)*subs(lambda1,t,x);
basisF(4) = subs(lambda1,t,y)*subs(lambda2,t,x);
basisF(5) = subs(lambda2,t,y)*subs(lambda2,t,x);
basisF(6) = subs(lambda3,t,y)*subs(lambda2,t,x);
basisF(7) = subs(lambda1,t,y)*subs(lambda3,t,x);
basisF(8) = subs(lambda2,t,y)*subs(lambda3,t,x);
basisF(9) = subs(lambda3,t,y)*subs(lambda3,t,x);

% localMass(i,j) = int_0^1 int_0^1 basisF(i)*basisF(j) dx dy
localMass = zeros(nrVBasisF,nrVBasisF); 
for i = 1:nrVBasisF
    for j = 1:nrVBasisF
        localMass(i,j) = int(int(basisF(i)*basisF(j),x,0,1),y,0,1);
    end
end

% calculate gradients
for i = 1:nrVBasisF
    basisFdx(i) = diff(basisF(i),x);
    basisFdy(i) = diff(basisF(i),y);
end

% for the diffusive term
% localStiff x_r x_t(i,j) = (int_0^1)^2 d/dx_r basisF(i)* d/dx_t basisF(j) dx dy
localStiff.xy = zeros(nrVBasisF,nrVBasisF);
localStiff.yx = zeros(nrVBasisF,nrVBasisF);
localStiff.xx = zeros(nrVBasisF,nrVBasisF);
localStiff.yy = zeros(nrVBasisF,nrVBasisF);
for i = 1:nrVBasisF
    for j = 1:nrVBasisF
        localStiff.xy(i,j) = int(int(basisFdx(i)*basisFdy(j),x,0,1),y,0,1);
        localStiff.yx(i,j) = int(int(basisFdy(i)*basisFdx(j),x,0,1),y,0,1);
		localStiff.xx(i,j) = int(int(basisFdx(i)*basisFdx(j),x,0,1),y,0,1);
		localStiff.yy(i,j) = int(int(basisFdy(i)*basisFdy(j),x,0,1),y,0,1);
    end
end

% nonlinear term (improve by using symmetries)
% localNonLin.x(k,r,l) = (int_0^1)^2 basisF(k) d/dx basisF(r) basisF(l) dx
localNonLin.x = zeros(nrVBasisF,nrVBasisF,nrVBasisF); 
localNonLin.y = zeros(nrVBasisF,nrVBasisF,nrVBasisF);
% for k = 1:nrVBasisF
%     for r = 1:nrVBasisF
%         for l = 1:nrVBasisF
%             localNonLin.x(k,r,l) = int(int(basisF(k)*basisFdx(r)*basisF(l),...
%                 x,0,1),y,0,1);
%             localNonLin.y(k,r,l) = int(int(basisF(k)*basisFdy(r)*basisF(l),...
%                 x,0,1),y,0,1);
%         end
%     end
% end

% pressure, localPdivV x_r(i,j) = int_0^1 int_0^1 basisP(j) d/dx_r basisF(i)
basisP(2) = x - 1/2;
basisP(3) = y - 1/2;
basisP(1) = 1;

localPdivV.x = zeros(nrPBasisF,nrVBasisF);
localPdivV.y = zeros(nrPBasisF,nrVBasisF);
for j = 1:nrVBasisF
    for i = 1:nrPBasisF
        localPdivV.x(i,j) = int(int(basisFdx(j)*basisP(i),x,0,1),y,0,1);
        localPdivV.y(i,j) = int(int(basisFdy(j)*basisP(i),x,0,1),y,0,1);
    end
end

localPMass = zeros(nrPBasisF, nrPBasisF);

for i = 1:nrPBasisF
    for j = 1:nrPBasisF
        localPMass(i,j) = int(int(basisP(j)*basisP(i),x,0,1),y,0,1);
    end
end
localMatrix = struct('vmass',localMass,'stiff', localStiff,...
     'pdivv',localPdivV, 'nonlin',localNonLin, 'pmass', localPMass);
end

