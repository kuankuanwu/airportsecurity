function [ snima,OCBA_solution,sencond_sema ] = OCBAf( sencond_tauavma,sencond_sema,deltaa )
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    global big_s1
    global big_s2
    global stage 
    
    nima = zeros(big_s1,big_s2);
    snima = zeros(big_s1,big_s2);
    cur_best_tau = 100;
    cur_best_s = 0;
    second_best_s=0;
    second_best_tau =100;
    OCBA_solution=[];


% find the best solution 
    for j=1:big_s1
        for k=1:big_s2

           if  sencond_tauavma(j+(k-1)*big_s1) < cur_best_tau
               cur_best_tau = sencond_tauavma(j+(k-1)*big_s1);
               cur_best_s = j+(k-1)*big_s1;
               best_s=[j,k];
            end
        end
    end

    if sencond_sema(cur_best_s)<=0
        sencond_sema(cur_best_s)=0.0001;
    end

%find the second best solution 
    for j=1:big_s1
        for k=1:big_s2

           if  j+(k-1)*big_s1 ~= cur_best_s  && sencond_sema(j+(k-1)*big_s1) > 0.0001
              if  sencond_tauavma(j+(k-1)*big_s1) < second_best_tau

                  second_best_tau = sencond_tauavma(j+(k-1)*big_s1);
                  second_best_s = j+(k-1)*big_s1;
                  second_s=[j,k];
              end
           end
        end
    end

    base = second_best_s;


    nima(base)=1;



    for j=1:big_s1
        for k=1:big_s2


            if (j+(k-1)*big_s1) ~= cur_best_s



                if sencond_tauavma(j+(k-1)*big_s1) < 1

                    if  (j+(k-1)*big_s1) ~= base

                       if sencond_sema(j+(k-1)*big_s1)<=0
                          sencond_sema(j+(k-1)*big_s1)=0.00001; 
                       end

                       if sencond_tauavma(j+(k-1)*big_s1) == sencond_tauavma(cur_best_s)
                           sencond_tauavma(j+(k-1)*big_s1) = sencond_tauavma(j+(k-1)*big_s1)+0.00001;
                       end    

                        nima(j+(k-1)*big_s1) = nima(base)*(sencond_sema(j+(k-1)*big_s1)*(sencond_tauavma(cur_best_s)-sencond_tauavma(base))/((sencond_sema(base))*(sencond_tauavma(cur_best_s)-sencond_tauavma(j+(k-1)*big_s1))))^2;


                    end
                end


            end
         end

    end





    nito=0;
    nitotal=0;

    for j=1:big_s1
        for k=1:big_s2
         if   j+(k-1)*big_s1 ~= cur_best_s && sencond_tauavma(j+(k-1)*big_s1) < 1

            nito =  nito +  (nima(j+(k-1)*big_s1)/sencond_sema(j+(k-1)*big_s1))^2;

         end
        end
    end

    nima(cur_best_s) =  sencond_sema(cur_best_s)*(nito)^0.5;
      

    for j=1:big_s1
        for k=1:big_s2
         if sencond_tauavma(j+(k-1)*big_s1) < 1
            nitotal= nitotal+nima(j+(k-1)*big_s1);
         end
        end
    end
         

    for j=1:big_s1
        for k=1:big_s2 
           if sencond_tauavma(j+(k-1)*big_s1) < 1
                snima(j+(k-1)*big_s1) = round(deltaa*nima(j+(k-1)*big_s1)/nitotal);
           end
           if snima(j+(k-1)*big_s1) > 0 
               solution=[stage,j,k,snima(j+(k-1)*big_s1),sencond_tauavma(j+(k-1)*big_s1),sencond_sema(j+(k-1)*big_s1)];
               OCBA_solution=[OCBA_solution;solution];
           end
        end
    end
    
%     if snima(cur_best_s)<=0
%         if sum(sum(snima)) >= 20
%             snima(cur_best_s) = 1;
%             solution=[stage,best_s,snima(cur_best_s)];
%             OCBA_solution=[OCBA_solution;solution];
%         else
%             snima(cur_best_s) = 20-sum(sum(snima));
%             solution=[stage,best_s,snima(cur_best_s)];
%             OCBA_solution=[OCBA_solution;solution];
%         end
%     else
%         if sum(sum(snima)) < 20
%             [a,~]=size(OCBA_solution);
%             for i = 1: a
%                 if sum(ismember(OCBA_solution(i,2:3),best_s)) == 2
%                     snima(cur_best_s) = snima(cur_best_s) + 20 - sum(sum(snima));
%                     OCBA_solution(i,4)= snima(cur_best_s);
%                 end
%             end
%            
%             
%         end
%     end
    
    if sum(sum(snima))==0
        snima(cur_best_s)=deltaa/2;
        snima(second_best_s)=deltaa/2;
        OCBA_solution=[stage,best_s,snima(cur_best_s);stage,second_s,snima(second_best_s)] ;
    end
end

