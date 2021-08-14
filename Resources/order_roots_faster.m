function ord_index = order_roots_faster(r,centers)
%DO NOT USE, NOT READY Order_roots: receives N roots and N centers,
%reorders the roots by considering centers one at a time
N = length(r);

ord_index=1:N;
for ii=1:N
    [~,indx]=min(abs(centers(ii)-r(ord_index(ii:N))));
    indx=indx+ii-1;
    tmp=ord_index(ii);
    ord_index(ii)=indx;
    ord_index(indx)=tmp;
end


end

