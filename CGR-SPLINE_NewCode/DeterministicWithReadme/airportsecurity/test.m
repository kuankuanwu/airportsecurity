iseedRet = zeros(1,RetailerNumber);
TmpOrderArrival_Ret = zeros(1,RetailerNumber);
OrderArrival_Ret = zeros(1,RetailerNumber);
Review_Ret = OrderArrival_War(1,1);
LastReview = Review_Ret ; 
point = 2;
MaxReviewTimes = 2;

for i = 1:RetailerNumber
iseedRet(i) = iseed(2)+1000*i;
end

for j = 1:RetailerNumber
    
[LeadTime,iseedRet(j)] =rexpo(mean,iseedRet(j));   
OrderArrival_Ret(1,j) =  OrderArrival_War(1,1)+ LeadTime;
end


 while (point <= MaxReviewTimes)

   if  OrderArrival_War(1,point) > LastReview && OrderArrival_War(1,point) < LastReview + tau
       Review_Ret = [Review_Ret;OrderArrival_War(1,point)];
       LastReview = OrderArrival_War(1,point) ; 
      
       for k = 1:RetailerNumber
           [LeadTime,iseedRet(k)] = rexpo(mean,iseedRet(k));
           TmpOrderArrival_Ret(1,k) = OrderArrival_War(1,point) + LeadTime;   
        end
       
       point = point + 1;
       
   else     
       Review_Ret = [Review_Ret ; LastReview+tau];
        
       for m = 1:RetailerNumber
           [LeadTime,iseedRet(m)] = rexpo(mean,iseedRet(m));
           TmpOrderArrival_Ret(1,m) = OrderArrival_War(1,point) + LeadTime;   
       end
        
       LastReview = LastReview+tau ;
   end
  
   OrderArrival_Ret = [ OrderArrival_Ret;  TmpOrderArrival_Ret] ;
   
 end