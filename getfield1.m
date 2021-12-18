function[field2]=getfield1(ph,voxelsize,padsize,mask)
% Authors: Janis Stiegeler, Sina Straub
% Email: j.stiegeler@dkfz.de, sina.straub@gmail.com, sina.straub@dkfz.de
% Date: 27.03.2021 V1.1

[uw, ~]=MRPhaseUnwrap(ph,'voxelsize',voxelsize,'padsize',padsize);
[tph,~]=V_SHARP(uw,mask,'voxelsize',voxelsize,'smvsize',12);

fieldin=((uw.*mask-tph.*mask));

[fieldx,fieldy,fieldz] = gradient((fieldin),voxelsize(1),voxelsize(2),voxelsize(3));

field = sqrt(fieldx.^2+fieldy.^2+fieldz.^2).*mask;
N=size(field);

se=strel('sphere',8);
maske=imerode(mask,se);
se=strel('sphere',12);

maske1 = imdilate(maske,se);

field_help = zeros(size(maske));
field_help(maske ==1) = field(maske ==1);

field1 = zeros(size(maske1));

for k = 1:N(1)
    for l = 1:N(2)
        for m = 1:N(3)
            if maske1(k,l,m)== 1
                for s = [3 6 12]
                    if (s<k) && (k <=N(1)-s) && (s<l) && (l<=N(2)-s) && (s<m) && (m<=N(3)-s)
                        index_help= field_help(k-s:k+s,l-s:l+s,m-s:m+s);
                        help_sum = sum(abs(index_help(:)));
                    end
                    
                    if (s >= 6) && (help_sum ~= 0)
                        break
                    end
                end
                index_help = index_help(index_help~=0);
                help_mean = mean(index_help(:));
                field1(k,l,m)= help_mean;
            end
            
            if (~(20<=k)) && (k<=N(1)-6)  && (7<=l) && (l<=N(2)-6) && (7<=m) && (m<=N(3)-6) %Rand 1a (Ränder ohne Ecken, könnte man noch machen, aber dort kein Feld)
                if maske1(k,l,m)== 1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k:k+6,l-6:l+6,m-6:m+6);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (k<=N(1)-a) && (a<l) && (l<=N(2)-a) && (a<m) && (m<=N(3)-a)
                                    index_help= field_help(k:k+a,l-a:l+a,m-a:m+a);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
            end
            if (7<=k) && (~(k<=N(1)-20)) && (7<=l) && (l<=N(2)-6) && (7<=m) && (m<=N(3)-6) %Rand 1b
                if maske1(k,l,m)==1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k-6:k,l-6:l+6,m-6:m+6);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (a<k) && (a<l) && (l<=N(2)-a) &&  (a<m) && (m<=N(3)-a)
                                    index_help= field_help(k-a:k,l-a:l+a,m-a:m+a);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
                
            end
            
            if (7<=k) && (k<=N(1)-6) && (~(20<=l)) && (l<=N(2)-6)&& (7<=m) && (m<=N(3)-6)  %Rand 2a
                if maske1(k,l,m)==1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k-6:k+6,l:l+6,m-6:m+6);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (a<k) && (k<=N(1)-a) && (l<=N(2)-a) &&  (a<m) && (m<=N(3)-a)
                                    index_help= field_help(k-a:k+a,l:l+a,m-a:m+a);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
            end
            if (7<=k) && (k<=N(1)-6) && (7<=l) && (~(l<=N(2)-20)) && (7<=m) && (m<=N(3)-6)  %Rans 2b
                if maske1(k,l,m)==1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k-6:k+6,l-6:l,m-6:m+6);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (a<k) && (k<=N(1)-a) && (a<l) &&  (a<m) && (m<=N(3)-a)
                                    index_help= field_help(k-a:k+a,l-a:l,m-a:m+a);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
            end
            if (7<=k) && (k<=N(1)-6) && (7<=l) && (l<=N(2)-6) && (~(20<=m)) && (m<=N(3)-6)  %Rans 3a
                if maske1(k,l,m)==1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k-6:k+6,l-6:l+6,m:m+6);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (a<k) && (k<=N(1)-a) && (a<l) && (l<=N(2)-a) &&  (m<=N(3)-a)
                                    index_help= field_help(k:k+a,l-a:l+a,m:m+a);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
            end
            if (7<=k) && (k<=N(1)-6) && (7<=l) && (l<=N(2)-6) && (7<=m) && (~(m<=N(3)-20))   %3b
                if maske1(k,l,m)==1
                    if (isnan(field1(k,l,m))) || (field1(k,l,m) == 0)
                        index_help= field_help(k-6:k+6,l-6:l+6,m-6:m);
                        if sum(abs(index_help(:))) == 0
                            for a = [7 8 10 12 14 16 18 20]
                                if (a<k) && (k<=N(1)-a) && (a<l) && (l<=N(2)-a) &&  (a<m)
                                    index_help= field_help(k:k+a,l-a:l+a,m-a:m);
                                end
                                if sum(abs(index_help(:))) ~= 0
                                    break
                                end
                            end
                        end
                    end
                    index_help = index_help(index_help~=0);
                    field1(k,l,m)=mean(index_help(:));
                end
            end
        end
    end
end

field = field1;

A= abs(field);
permax = prctile(A(:),99.9);
permin = prctile(A(:),90);
field2 = zeros(size(field));
field2(A<=permax) = ((permax-A(A<=permax))/(permax-permin));
field2(A<permin) = 1;
field2 = field2.*mask;

