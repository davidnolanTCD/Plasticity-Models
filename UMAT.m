function [ddsdde, stress, statev] = UMAT(stress, dstran, stran, props, statev)
disp('--------------UMAT START--------------')
%     Specify material properties

      E = props(1);
      xnue = props(2);
      sigy0 = props(3);
      h = props(4);

%    Recover effective plastic strain, p, and isotropic 
%    hardening variable, r,from previous time step
      p = statev(1);
      r = statev(2);
      toler = 10^-4;
      
%    Properties

      ebulk3 = E/(1-2*xnue);
      xk = ebulk3/3;
      eg2 = E/(1+xnue);
      eg = eg2/2;


%    Set up 6x6 elasticity matrix

       [ddsdde] = Celastic(E, xnue);

%     Save stress at beginning of time step in stressold

      stressold = stress;
      
%     Obtain trial (elastic) stress

      dstress = ddsdde*dstran;      
      stress = stress + dstress;
      
%     Write trial stresses in matrix form
      [str] = tensizer(stress);
      

%     Calculate deviatoric trial stress
 
      devstr = str - (trace(str)*eye(3,3)./3);

%     Calculate effective trial stress

      [pj] = vonmises(stress);
      
%    determine trial flow direction for use with the jacobian

      xndir = devstr./pj;

%    ...and write in voigt notation

        [xnv] = vectorizer(xndir);
        

%    determine if the yield condition is satisfied

      zy = pj - r - sigy0;

      if (zy > 0.0)
          
          disp('plasticity')

%     Use newton iteration to determine effective plastic strain increment
 
          r0 = r;
          dp = 0.0;
      
              for n = 1:10
                  
                  res = pj-3.*eg*dp-r-sigy0;
                  dp = dp + res/(3.*eg+h);
                  r = r0 + h*dp;
              
                  if(abs(res) < toler)
                    break
                  end
              end
              
              if(n==10)
                 disp('warning max iters exceeded') 
              end

%     Determine the increments in plastic strain


       dpstrn = (3/2)*dp*(devstr./pj);

%     Write the strain increments in voigt notation (with engg shears)

        dpstran = vectorizer(dpstrn);
        
        for i = 4:6
            dpstran(i) = 2*dpstran(i);
        end
        
      
%     Calculate the elastic strain increments            
      
      destran = dstran - dpstran;

%     Determine stress increment

        dstress = ddsdde*destran;

%     Update the stress, effective plastic strain
%     (note: isotropic hardening variable already updated)

       stress = stressold + dstress;

       p = p + dp;
       
%     store updated state variables

       statev(1) = p;
       statev(2) = r;

%%    Determine jacobian

     ddsdde = zeros(6,6);
     xiden = eye(3,3);
     

      xr = (pj-3.*eg*dp)/pj;
      q = (1/(1+3.*eg/h)-xr)*3/2;
      for i=1:3
        for j=1:3
            ddsdde(i,j) = 2*eg*q*xnv(i)*xnv(j)...
            + (xk-eg*xr*2/3) + 2*eg*xr*xiden(i,j);
        end 
      end
      
      for l = 4:6
          
        for k = 1:3
            ddsdde(k,l) = 2*eg*q*xnv(k)*xnv(l);
            ddsdde(l,k) = ddsdde(k,l);
        end
        
        ddsdde(l,l) = 2*eg*q*xnv(l)*xnv(l) + eg*xr;
      end
      
      
      end % End of if statement
disp('-------------UMAT END--------------')
end