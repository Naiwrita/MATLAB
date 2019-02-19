%Start Programming.
clc
clear all
close all

%Resources...
foodr = input('Enter available food(Ton): \n');
waterr = input('Enter available water(Gallon): \n');
clothr = input('Enter available cloth(Dozen): \n');
mer = input('Enter available medical equipments(Packet): \n');
rmr = input('Enter available rescue member(Number): \n');
s=input('How many sites? ');


totfoodd=0;
totwaterd=0;
totclothd=0;
totmed=0;
totrmd=0;


%initial details of each area...
for n = 1:s
   %Initial demand
    ifdmnd(n) = input(['Enter food Demand of site ' num2str(n) '(Ton):\n']);
    iwdmnd(n) = input(['Enter water Demand of site ' num2str(n) '(Gallon):\n']);
    icdmnd(n) = input(['Enter cloths Demand of site ' num2str(n) '(Dozen):\n']);
    imedmnd(n) = input(['Enter medical equipments Demand of site ' num2str(n) '(Packet):\n']);
    irmdmnd(n) = input(['Enter rescue member Demand of site ' num2str(n) '(Number):\n']);
        
    %calculate the total resource demand for each site
    totfoodd=totfoodd+ifdmnd(n);
    totwaterd=totwaterd+iwdmnd(n);
    totclothd=totclothd+icdmnd(n);
    totmed=totmed+imedmnd(n);
    totrmd=totrmd+irmdmnd(n);
end
minfdmnd = min(ifdmnd);
minwdmnd = min(iwdmnd);
mincdmnd = min(icdmnd);
minmdmnd = min(imedmnd);
minrdmnd = min(irmdmnd);

 totpdf = 0;
 totpdw = 0;
 totpdc = 0;
 totpdm = 0;
 totpdr = 0;
 
%finding plotted demand
 for i=1:s
    fdmnd(i)=(ifdmnd(i)/totfoodd)*foodr; 
    totpdf = totpdf+fdmnd(i);
    wdmnd(i)=(iwdmnd(i)/totwaterd)*waterr;
    totpdw = totpdw+wdmnd(i);
    cdmnd(i)=(icdmnd(i)/totclothd)*clothr;
    totpdc = totpdc+cdmnd(i);
    medmnd(i)=(imedmnd(i)/totmed)*mer;
    totpdm = totpdm+medmnd(i);
    rmdmnd(i)=(irmdmnd(i)/totrmd)*rmr;
    totpdr = totpdr+rmdmnd(i);
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



k = 1;
if totwaterd >= waterr
  for m = 1:p1+1
    for n =1:s
       IWCS(m,n) = round(1/(m+1).*wdmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        IWCS(i,j) = round(wdmnd(j)-(k*5));
    end
    k = k+1;
end
else
  for m = 1:p1+1
    for n =1:s
        IWCS(m,n) = round(1/(m+1).*iwdmnd(n));
    end
    m = m+1;
  end
  for i = p1:p
    for j = 1:s
        IWCS(i,j) = round(iwdmnd(j)-(k*5));
    end
    k = k+1;
    end
end
k = 1;
if totclothd >= clothr
  for m = 1:p1+1
    for n =1:s
        ICCS(m,n) = round(1/(m+1).*cdmnd(n));
    end
    m = m+1;
  end
   for i = p1+1:p
    for j = 1:s
        ICCS(i,j) = round(cdmnd(j)-(k*5));
    end
    k = k+1;
end
else
  for m = 1:p1+1
    for n =1:s
        ICCS(m,n) = round(1/(m+1).*icdmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        ICCS(i,j) = round(icdmnd(j)-(k*5));
    end
    k = k+1;
end
end
 
k = 1;
if totmed >= mer
  for m = 1:p1+1
    for n =1:s
        IMCS(m,n) = round(1/(m+1).*medmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        IMCS(i,j) = round(medmnd(j)-(k*5));
    end
    k = k+1;
end
else
  for m = 1:p1+1
    for n =1:s
         IMCS(m,n) = round(1/(m+1).*imedmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        IMCS(i,j) = round(imedmnd(j)-(k*5));
    end
    k = k+1;
end
end

k = 1;
if totrmd >= rmr
  for m = 1:p1+1
    for n =1:s
        IRCS(m,n) = round(1/(m+1).*rmdmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        IRCS(i,j) = round(rmdmnd(j)-(k*5));
    end
    k = k+1;
end
else
  for m = 1:p1+1
    for n =1:s
         IRCS(m,n) = round(1/(m+1).*irmdmnd(n));
    end
    m = m+1;
  end
  for i = p1+1:p
    for j = 1:s
        IRCS(i,j) = round(irmdmnd(j)-(k*5));
    end
    k = k+1;
end
end


%Displaying the details of initial candidate solutions 
disp('Initial food allocations');
disp(IFCS);
disp('Initial water allocations');
disp(IWCS);
disp('Initial cloth allocations');
disp(ICCS);
disp('Initial medicine allocations');
disp(IMCS);
disp('Initial rescueteam allocations'); 
disp(IRCS);


 kf = 1;
 kw = 1;
 kc = 1;
 km = 1;
 kr = 1;
%{
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
    
    
    for i = 1:p
         totwcs = 0;
        for j = 1:s
             totwcs = totwcs+IWCS(i,j);
             j = j+1;
        end
        if totwaterd < waterr
            fitwcs = totwaterd-totwcs;
        else
            fitwcs = totpdw-totwcs;
        end
        IWCS(i,s+1) = fitwcs;
        i=i+1;
     end
    tot = 0;
    tot = sum(IWCS(:,s+1));
    FITWT(kw) = tot/p;
    kw = kw+1;
    
    
    for i = 1:p
           totccs = 0;
        for j = 1:s
           totccs = totccs+ICCS(i,j);
           j = j+1;
        end
        if totclothd < clothr
            fitccs = totclothd-totccs;
        else
            fitccs = totpdc-totccs;
        end
        ICCS(i,s+1) = fitccs;
        i=i+1;
      end
      tot = 0;
      tot = sum(ICCS(:,s+1));
      FITFD(kc) = tot/p;
      kc = kc+1;
      
       
       for i = 1:p
            totmcs = 0;
            for j = 1:s
                totmcs = totmcs+IMCS(i,j);
                j = j+1;
            end
            if totmed < mer
               fitmcs = totmed-totmcs;
            else
               fitmcs = totpdm-totmcs;
            end
            IMCS(i,s+1) = fitmcs;
            i=i+1;
       end
       tot = 0;
       tot = sum(IMCS(:,s+1));
       FITMD(km) = tot/p;
       km = km+1;
    
       
        for i = 1:p
             totrcs = 0;
            for j = 1:s
                totrcs = totrcs+IRCS(i,j);
                j = j+1;
            end
             if totrmd <= rmr
                fitrcs = totrmd-totmcs;
             else
                fitrcs = totpdr-totmcs;
             end
            IRCS(i,s+1) = fitrcs;
            i=i+1;
        end
        tot = 0;
        tot = sum(IRCS(:,s+1));
        FITRS(kr) = tot/p;
        kr = kr+1;
           
    [fdmax,fdmaxpos] = max(IFCS(:,s+1));
    [wtmax,wtmaxpos] = max(IWCS(:,s+1));
    [clmax,clmaxpos] = max(ICCS(:,s+1));
    [medmax,medmaxpos] = max(IMCS(:,s+1));
    [rsmax,rsmaxpos] = max(IRCS(:,s+1));
    
    [fdmin,fdminpos] = min(IFCS(:,s+1));
    [wtmin,wtminpos] = min(IWCS(:,s+1));
    [clmin,clminpos] = min(ICCS(:,s+1));
    [medmin,medminpos] = min(IMCS(:,s+1));
    [rsmin,rsminpos] = min(IRCS(:,s+1));
  
      %sorting population matrix based on fitness value
      IFCSS = sortrows(IFCS, (s+1));
      IWCSS = sortrows(IWCS, (s+1));
      ICCSS = sortrows(ICCS, (s+1));
      IMCSS = sortrows(IMCS, (s+1));
      IRCSS = sortrows(IRCS, (s+1));
        
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
     
      %WATER
      MPIWCS = IWCS;
      for j = 1:s+1
        MPIWCS(wtmaxpos,j) = IWCS(wtminpos,j);
      end
      MPIWCS = sortrows(MPIWCS, (s+1));
     
      for i = 1:p1
          rndw(i) = -1;
          i = i+1;
      end
      i=1;
      while i<=p1
          flag=0;
          r = randi([p1+1,p],1,1);
          for j= 1:i
              if(r==rndw(j))
                  flag = 1;
                  break;
              end
              j=j+1;
          end
          if (flag~=1)
              rndw(i)=r;
              i=i+1;
          end
      end
      
      for i = 1:p1
          MPIWCS(i,s+2) = rndw(i);
          W=rndw(i);
          MPIWCS(W,s+2) = i;
          i = i+1;
      end
      
      
      %CLOTH
      MPICCS = ICCS;
      for j = 1:s+1
        MPICCS(clmaxpos,j) = ICCS(clminpos,j);
      end
      MPICCS = sortrows(MPICCS, (s+1));
      
      for i = 1:p1
          rndc(i) = -1;
          i = i+1;
      end
      i=1;
      while i<=p1
          flag=0;
          r = randi([p1+1,p],1,1);
          for j= 1:i
              if(r==rndc(j))
                  flag = 1;
                  break;
              end
              j=j+1;
          end
          if (flag~=1)
              rndc(i)=r;
              i=i+1;
          end
      end
      
      for i = 1:p1
          MPICCS(i,s+2) = rndc(i);
          W=rndc(i);
          MPICCS(W,s+2) = i;
          i = i+1;
      end
     
      %MEDICAL
      MPIMCS = IMCS;
      for j = 1:s+1
        MPIMCS(medmaxpos,j) = IMCS(medminpos,j);
      end
      MPIMCS = sortrows(MPIMCS, (s+1));
      
      for i = 1:p1
          rndm(i) = -1;
          i = i+1;
      end
      i=1;
      while i<=p1
          flag=0;
          r = randi([p1+1,p],1,1);
          for j= 1:i
              if(r==rndm(j))
                  flag = 1;
                  break;
              end
              j=j+1;
          end
          if (flag~=1)
              rndm(i)=r;
              i=i+1;
          end
      end
      
      for i = 1:p1
          MPIMCS(i,s+2) = rndm(i);
          W=rndm(i);
          MPIMCS(W,s+2) = i;
          i = i+1;
      end
     
      
      %RESCUE
      MPIRCS = IRCS;
      for j = 1:s+1
        MPIRCS(rsmaxpos,j) = IRCS(rsminpos,j);
      end
      MPIRCS = sortrows(MPIRCS, (s+1));
          i = i+1;
      end
      i=1;
      while i<=p1
          flag=0;
          r = randi([p1+1,p],1,1);
          for j= 1:i
      
      for i = 1:p1
          rndr(i) = -1;
              if(r==rndr(j))
                  flag = 1;
                  break;
              end
              j=j+1;
          end
          if (flag~=1)
              rndr(i)=r;
              i=i+1;
          end
      end
      
      for i = 1:p1
          MPIRCS(i,s+2) = rndr(i);
          W=rndr(i);
          MPIRCS(W,s+2) = i;
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
  
      if(C_Rate == 1 || C_Rate == 0.9) %all parents undergo crossover
          for  i = 1:p1
              for j = 1:s
                  b = MPIWCS(i,s+2);
                  OFSPRING1(i,j) = a.*MPIWCS(i,j)+(1-a).*MPIWCS(b,j);
                  OFSPRING2(i,j) = (1-a).*MPIWCS(i,j)+a.*MPIWCS(b,j);
              end
          end
      else
          for  i = 1:p1  
              for j = 1:s
                    b = MPIWCS(p,s+2);
                    OFSPRING1(i,j) = round(a.*MPIWCS(p,j)+(1-a).*MPIWCS(b,j));
                    OFSPRING2(i,j) =round((1-a).*MPIWCS(p,j)+a.*MPIWCS(b,j));
              end
         end
          
      end
      
      MPIWCS=[OFSPRING1;OFSPRING2];
      
      if(C_Rate == 1 || C_Rate == 0.9) %all parents undergo crossover
          for  i = 1:p1
              for j = 1:s
                  b = MPICCS(i,s+2);
                  OFSPRING1(i,j) = a.*MPICCS(i,j)+(1-a).*MPICCS(b,j);
                  OFSPRING2(i,j) = (1-a).*MPICCS(i,j)+a.*MPICCS(b,j);
              end
          end
      else
          for  i = 1:p1  
              for j = 1:s
                    b = MPICCS(p,s+2);
                    OFSPRING1(i,j) = round(a.*MPICCS(p,j)+(1-a).*MPICCS(b,j));
                    OFSPRING2(i,j) =round((1-a).*MPICCS(p,j)+a.*MPICCS(b,j));
              end
         end
      end
      MPICCS=[OFSPRING1;OFSPRING2];
      
      if(C_Rate == 1 || C_Rate == 0.9) %all parents undergo crossover
          for  i = 1:p1
              for j = 1:s
                  b = MPIMCS(i,s+2);
                  OFSPRING1(i,j) = a.*MPIMCS(i,j)+(1-a).*MPIMCS(b,j);
                  OFSPRING2(i,j) = (1-a).*MPIMCS(i,j)+a.*MPIMCS(b,j);
              end
          end
      else
          for  i = 1:p1  
              for j = 1:s
                    b = MPIMCS(p,s+2);
                    OFSPRING1(i,j) = round(a.*MPIMCS(p,j)+(1-a).*MPIMCS(b,j));
                    OFSPRING2(i,j) =round((1-a).*MPIMCS(p,j)+a.*MPIMCS(b,j));
              end
         end
      end
      MPIMCS=[OFSPRING1;OFSPRING2];
      
      if(C_Rate == 1 || C_Rate == 0.9) %all parents undergo crossover
          for  i = 1:p1
              for j = 1:s
                  b = MPIRCS(i,s+2);
                  OFSPRING1(i,j) = a.*MPIRCS(i,j)+(1-a).*MPIRCS(b,j);
                  OFSPRING2(i,j) = (1-a).*MPIRCS(i,j)+a.*MPIRCS(b,j);
              end
          end
      else
          for  i = 1:p1  
              for j = 1:s
                    b = MPIRCS(p,s+2);
                    OFSPRING1(i,j) = round(a.*MPIRCS(p,j)+(1-a).*MPIRCS(b,j));
                    OFSPRING2(i,j) =round((1-a).*MPIRCS(p,j)+a.*MPIRCS(b,j));
              end
         end
      end
      MPIRCS=[OFSPRING1;OFSPRING2];
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
              t= MPIWCS(i,pos1);
              MPIWCS(i,pos1)=MPIWCS(i,pos2);
              MPIWCS(i,pos2)=t;
              c=c+1;
          end
          if(c==PM)
              break;
          end
      end
      
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
              t= MPICCS(i,pos1);
              MPICCS(i,pos1)=MPICCS(i,pos2);
              MPICCS(i,pos2)=t;
              c=c+1;
          end
          if(c==PM)
              break;
          end
      end
      
      
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
              t= MPIMCS(i,pos1);
              MPIMCS(i,pos1)=MPIMCS(i,pos2);
              MPIMCS(i,pos2)=t;
              c=c+1;
          end
          if(c==PM)
              break;
          end
      end
      
      
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
              t= MPIRCS(i,pos1);
              MPIRCS(i,pos1)=MPIRCS(i,pos2);
              MPIRCS(i,pos2)=t;
              c=c+1;
          end
          if(c==PM)
              break;
          end
      end
     
      IFCS=MPIFCS;
      IWCS=MPIWCS;
      ICCS=MPICCS;
      IMCS=MPIMCS;
      IRCS=MPIRCS;
      
       if(S==150)
           break;
       end
      S=S+1;
 end
       
   
 %}         
         
          
     

                
                            
              
      
      
      
      
          
          
      
      
          
          
          
          
          
         
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
      
      
      
      