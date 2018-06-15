clear all


A=34;      %degree

X=123;

Y=-32;

Z=23;

normalization_factor=(X^2+Y^2+Z^2)^0.5;

vector=[X Y Z]/normalization_factor;

X=vector(1);
Y=vector(2);
Z=vector(3);

theta=A/180*pi;


R1=[cos(theta)+(1-cos(theta))*X^2   -Z*sin(theta)+X*Y*(1-cos(theta))   Y*sin(theta)+X*Z*(1-cos(theta));  Z*sin(theta)+X*Y*(1-cos(theta))   cos(theta)+(1-cos(theta))*Y^2   -X*sin(theta)+Y*Z*(1-cos(theta));  -Y*sin(theta)+X*Z*(1-cos(theta))   X*sin(theta)+Y*Z*(1-cos(theta))   cos(theta)+(1-cos(theta))*Z^2];


R1(abs(R1)<0.00000000000001)=0;


vector=[0 0 1];
X=vector(1);
Y=vector(2);
Z=vector(3);
theta=90/180*pi;

R2=[cos(theta)+(1-cos(theta))*X^2   -Z*sin(theta)+X*Y*(1-cos(theta))   Y*sin(theta)+X*Z*(1-cos(theta));  Z*sin(theta)+X*Y*(1-cos(theta))   cos(theta)+(1-cos(theta))*Y^2   -X*sin(theta)+Y*Z*(1-cos(theta));  -Y*sin(theta)+X*Z*(1-cos(theta))   X*sin(theta)+Y*Z*(1-cos(theta))   cos(theta)+(1-cos(theta))*Z^2];


R2(abs(R2)<0.00000000000001)=0;

Rtotal=R2*R1;

Rtotal(abs(Rtotal)<0.00000000000001)=0;

Zsintheta2=Rtotal(2,1)-Rtotal(1,2);

Ysintheta2=Rtotal(1,3)-Rtotal(3,1);

Xsintheta2=Rtotal(3,2)-Rtotal(2,3);

temp1=1000000;
temp2=1000000;

for Z=-1:0.001:1
    for theta=0:pi/10000:pi
        A=Zsintheta2-2*Z*sin(theta);
        B=Rtotal(3,3)-(cos(theta)+(1-cos(theta))*Z^2);
        if ((A^2)+(B^2)<temp1) && (theta <= pi/2)
            temp1=(A^2)+(B^2);
            temp_theta1=theta;
            temp_Z1=Z;
        end
        if ((A^2)+(B^2)<temp2) && (theta > pi/2)
            temp2=(A^2)+(B^2);
            temp_theta2=theta;
            temp_Z2=Z;
        end
    end
end

Z_req1=temp_Z1;
Z_req2=temp_Z2;

theta_req1=temp_theta1;
theta_req2=temp_theta2;

temp1=1000000;
temp2=1000000;

for X=-1:0.0001:1
        A1=Xsintheta2-2*X*sin(theta_req1);
        B1=Rtotal(1,1)-(cos(theta_req1)+(1-cos(theta_req1))*X^2);
        if (A1^2)+(B1^2)<temp1
            temp1=(A1^2)+(B1^2);
            temp_X1=X;
        end
        
        A2=Xsintheta2-2*X*sin(theta_req2);
        B2=Rtotal(1,1)-(cos(theta_req2)+(1-cos(theta_req2))*X^2);
        if (A2^2)+(B2^2)<temp2
            temp2=(A2^2)+(B2^2);
            temp_X2=X;
        end
end

X_req1=temp_X1;
X_req2=temp_X2;

temp1=1000000;
temp2=1000000;

for Y=-1:0.0001:1
        A1=Ysintheta2-2*Y*sin(theta_req1);
        B1=Rtotal(2,2)-(cos(theta_req1)+(1-cos(theta_req1))*Y^2);
        if (A1^2)+(B1^2)<temp1
            temp1=(A1^2)+(B1^2);
            temp_Y1=Y;
        end
        A2=Ysintheta2-2*Y*sin(theta_req2);
        B2=Rtotal(2,2)-(cos(theta_req2)+(1-cos(theta_req2))*Y^2);
        if (A2^2)+(B2^2)<temp2
            temp2=(A2^2)+(B2^2);
            temp_Y2=Y;
        end
end

Y_req1=temp_Y1;
Y_req2=temp_Y2;




X=X_req1;
Y=Y_req1;
Z=Z_req1;
theta=theta_req1;

R_req1=[cos(theta)+(1-cos(theta))*X^2   -Z*sin(theta)+X*Y*(1-cos(theta))   Y*sin(theta)+X*Z*(1-cos(theta));  Z*sin(theta)+X*Y*(1-cos(theta))   cos(theta)+(1-cos(theta))*Y^2   -X*sin(theta)+Y*Z*(1-cos(theta));  -Y*sin(theta)+X*Z*(1-cos(theta))   X*sin(theta)+Y*Z*(1-cos(theta))   cos(theta)+(1-cos(theta))*Z^2];

residual_1=R_req1-Rtotal;



X=X_req2;
Y=Y_req2;
Z=Z_req2;
theta=theta_req2;

R_req2=[cos(theta)+(1-cos(theta))*X^2   -Z*sin(theta)+X*Y*(1-cos(theta))   Y*sin(theta)+X*Z*(1-cos(theta));  Z*sin(theta)+X*Y*(1-cos(theta))   cos(theta)+(1-cos(theta))*Y^2   -X*sin(theta)+Y*Z*(1-cos(theta));  -Y*sin(theta)+X*Z*(1-cos(theta))   X*sin(theta)+Y*Z*(1-cos(theta))   cos(theta)+(1-cos(theta))*Z^2];

residual_2=R_req2-Rtotal;