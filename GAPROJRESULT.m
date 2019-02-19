foodr = input('Enter available food(Ton): \n');
s=input('How many sites? ');

totfoodd=0;

%initial details of each area...
for n = 1:s
   %Initial demand
    ifdmnd(n) = input(['Enter food Demand of site ' num2str(n) '(Ton):\n']);
    %calculate the total resource demand for each site
    totfoodd=totfoodd+ifdmnd(n);
end
minfdmnd = min(ifdmnd);
totpdf = 0;
 
 
%finding plotted demand [this limits the allocations and does not allow the total allocation solution to exceed the availability]
 for i=1:s
    fdmnd(i)=(ifdmnd(i)/totfoodd)*foodr; 
    totpdf = totpdf+fdmnd(i);
 end
 
%forming initial candidate solutions
p = 2*s;
p1 = p/2;
k=1;
if totfoodd >= foodr
  for m = 1:p-1
    for n =1:s
        IFCS(m,n) = round(1/m.*fdmnd(n));
    end
    m = m+1;
  end
  for i = p-1:p
    for j = 1:s
        IFCS(i,j) = minfdmnd;
        
    end
    k = k+1;
end
else
  for m = 1:p-1
    for n =1:s
        IFCS(m,n) = round(1/m.*ifdmnd(n));
    end
    m = m+1;
  end
   for i = p-1:p
    for j = 1:s
        IFCS(i,j)= minfdmnd;
    end
    k = k+1;
    end
end

%Displaying the details of initial candidate solutions 
disp('Initial food allocations');
disp(IFCS);
kf = 1;
S=1;
 while 1
   %finding total allocations for each resource by each candidate solution
   tot = 0;
   for i = 1:p
         totfcs = 0;
        for j = 1:s
            totfcs = totfcs+IFCS(i,j);
        end
        if totfoodd < foodr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fitfcs = totfoodd-totfcs;
        else
            fitfcs = totpdf-totfcs;
        end
        IFCS(i,s+1) = fitfcs;
        tot = tot+fitfcs;
        i=i+1;
    end
    FITFD(k) = tot/p;
    kf = kf+1;
    
    
  
           
    [fdmax,fdmaxpos] = max(IFCS(:,s+1));
    [fdmin,fdminpos] = min(IFCS(:,s+1));
    
  
      %sorting population matrix based on fitness value
      IFCSS = sortrows(IFCS, (s+1));
    
        
      %MATING POOL AND CROSSOVER PARTNER
      %FOOD
      MPIFCS = IFCS;
      for j = 1:s+1
        MPIFCS(fdmaxpos,j) = IFCS(fdminpos,j);
      end
      MPIFCS = sortrows(MPIFCS, (s+1)); %in ascending order
   
      p1=p/2;
      for i = 1:p1
          rndf(i) = -1;
          i = i+1;
      end
      i=1;
      while i<=p1
          flag=0;
          r = randi([p1+1,p],1,1);
          for j= 1:i
              if(r==rndf(j))
                  flag = 1;
                  break;
              end
              j=j+1;
          end
          if (flag~=1)
              rndf(i)=r;
              i=i+1;
          end
      end
      
      for i = 1:p1
          MPIFCS(i,s+2) = rndf(i);
          W=rndf(i);
          MPIFCS(W,s+2) = i;
          i = i+1;
      end
     
      
     
      %CROSSOVER BASED ON CROSSOVER RATE
      a = 0.5;
      C_Rate = 0.8;
      PC = C_Rate*p;
      if(C_Rate == 1 || C_Rate == 0.9) %all parents undergo crossover
          for  i = 1:p1
              for j = 1:s
                  b = MPIFCS(i,s+2);
                  OFSPRING1(i,j) = round(a.*MPIFCS(i,j)+(1-a).*MPIFCS(b,j));
                  OFSPRING2(i,j) = round((1-a).*MPIFCS(i,j)+a.*MPIFCS(b,j));
              end
          end
      else
         for  i = 1:p1  
              for j = 1:s
                    b = MPIFCS(p,s+2);
                    OFSPRING1(i,j) = round(a.*MPIFCS(p,j)+(1-a).*MPIFCS(b,j));
                    OFSPRING2(i,j) =round((1-a).*MPIFCS(p,j)+a.*MPIFCS(b,j));
              end
         end
      end
      MPIFCS=[OFSPRING1;OFSPRING2];    
  
      
      %%%%%%%%%%%%%%%%%%%% try out other crossover rates
  
      %Mutation
     
      M_Rate = 0.01;
      PM = ceil(M_Rate*p*s);
      
      for i = 1:p
          X=0+(1-0).*rand(1,1);
          c = 0;
          if(X<=0.05)
              pos1 = randi(s);
              if(pos1 == s)
                    pos2 = 1;
              else
                    pos2 = pos1+1;
              end
              t= MPIFCS(i,pos1);
              MPIFCS(i,pos1)=MPIFCS(i,pos2);
              MPIFCS(i,pos2)=t;
              c=c+1;
          end
          if(c==PM)
              break;
          end
      end
    
      
     
      IFCS=MPIFCS;
     
       if(S==150)
           break;
       end
      S=S+1;
 end
       
          
         
          
     

                
                            
              
      
      
      
      
          
          
      
      
          
          
          
          
          
         
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
      
      
      
      
