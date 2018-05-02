%--------------------------------------------------------------------------
%            Lossless JPEG Compression 
%--------------------------------------------------------------------------


function compress_image

orig_im=imread('/home/utkumukan/Masaüstü/resim/1.jpg');

tmp_size=size(orig_im);


comp_R = orig_im(:,:,1);
comp_G = orig_im(:,:,2);
comp_B = orig_im(:,:,3);

ARed=padarray(comp_R,[1,0],0);
AGreen=padarray(comp_G,[1,0],0);
ABlue=padarray(comp_B,[1,0],0);

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Red Channel
for i = 2:size(ARed,1)-1
    for j = 2:size(ARed,2)
        fx=round(ARed(i,j-1)+((ARed(i-1,j)-ARed(i-1,j-1)))/2);              % P5: A + (B - C) / 2 
         comp_R(i-1,j)=comp_R(i-1,j)-fx;
    end
end

%Green Channel
for i = 2:size(AGreen,1)-1
    for j = 2:size(AGreen,2)
        fx=round(AGreen(i,j-1)+((AGreen(i-1,j)-AGreen(i-1,j-1)))/2);        % P5: A + (B - C) / 2 
         comp_G(i-1,j)=comp_G(i-1,j)-fx;
    end
end

%Blue Channel
for i = 2:size(ABlue,1)-1
    for j = 2:size(ABlue,2)
        fx=round(ABlue(i,j-1)+((ABlue(i-1,j)-ABlue(i-1,j-1)))/2);           % P5: A + (B - C) / 2 
         comp_B(i-1,j)=comp_B(i-1,j)-fx;
    end
end
%--------------------------------------------------------------------------
%Huffman Tree
%--------------------------------------------------------------------------
treeR=build_huffman(imhist(comp_R));
treeG=build_huffman(imhist(comp_G));
treeB=build_huffman(imhist(comp_B));
disp('->HUFFMAN');
%--------------------------------------------------------------------------
%Traversing
%--------------------------------------------------------------------------
sR=traverse(treeR);
sG=traverse(treeG);
sB=traverse(treeB);
disp('->TRAVERSE');

streamR=encode(comp_R,sR);
streamG=encode(comp_G,sG);
streamB=encode(comp_B,sB);
disp('->ENCODE');


retrieved_R=decode(streamR,treeR,tmp_size(1),tmp_size(2));
retrieved_G=decode(streamG,treeG,tmp_size(1),tmp_size(2));
retrieved_B=decode(streamB,treeB,tmp_size(1),tmp_size(2));
disp('->DECODE');

%--------------------------------------------------------------------------
%Decompression operation
%-------------------------------------------------------------------------
retrieved_R=uint8(retrieved_R);
retrieved_G=uint8(retrieved_G);
retrieved_B=uint8(retrieved_B);


AARed=padarray(retrieved_R,[1,0],0);
AAGreen=padarray(retrieved_G,[1,0],0);
AABlue=padarray(retrieved_B,[1,0],0);

%Red Channel
for i = 2:size(AARed,1)-1
    for j = 2:size(AARed,2)
        fx=round(ARed(i,j-1)+((AARed(i-1,j)-AARed(i-1,j-1)))/2);              % P5: A + (B -C) / 2
        retrieved_R(i-1,j)=retrieved_R(i-1,j)+fx;
        AARed=padarray(retrieved_R,[1,0],0);
    end
end

%Green Channel
for i = 2:size(AAGreen,1)-1
    for j = 2:size(AAGreen,2)
        fx=round(AGreen(i,j-1)+((AAGreen(i-1,j)-AAGreen(i-1,j-1)))/2);        % P5: A + (B - C) / 2
        retrieved_G(i-1,j)=retrieved_G(i-1,j)+fx;
        AAGreen=padarray(retrieved_G,[1,0],0);
    end
end

%Blue Channel
for i = 2:size(AABlue,1)-1
    for j = 2:size(AABlue,2)
        fx=round(ABlue(i,j-1)+((AABlue(i-1,j)-AABlue(i-1,j-1)))/2);           % P5: A + (B � C) / 2
        retrieved_B(i-1,j)=retrieved_B(i-1,j)+fx;
        AABlue=padarray(retrieved_B,[1,0],0);
    end
end


retrieved(tmp_size(1),tmp_size(2),tmp_size(3))=0;
retrieved(:,:,1)=retrieved_R;
retrieved(:,:,2)=retrieved_G;
retrieved(:,:,3)=retrieved_B;
retrieved=uint8(retrieved);

conc = (tmp_size(1)*tmp_size(2)*24)/(length(streamR)+length(streamG)+length(streamB));

figure
imshow(orig_im);
title('Original Image');
figure
imshow(retrieved);
title({'Compressed Image';'Compression Ratio :'; conc} );

imwrite(retrieved,'/home/utkumukan/Masaüstü/resim/1.jpg');



disp('------------------------------------------------------------------');
fprintf('Compression Ratio = %f\n',(tmp_size(1)*tmp_size(2)*24)/(length(streamR)+length(streamG)+length(streamB)));
disp('------------------------------------------------------------------');
disp('FINISHED');

end

%--------------------------------------------------------------------------
%Decoding Function
%-------------------------------------------------------------------------
function ret=decode(stream,tree,x,y)

ret(x,y)=0;
top=0;
for i=1:length(tree)
    if tree(i).prnt==-1
        top=i;
        break;
    end
end

i=1;
r=1;
c=1;
j=top;

while i<=length(stream)
    
    if (tree(j).char~=-1)
        ret(r,c)=tree(j).char-1;
        c=c+1;
        if c>y
            c=1;
            r=r+1;
        end
        j=top;
    end
    
    if stream(i)==1
        j=tree(j).c1;
    else
        j=tree(j).c2;
    end
    
    i=i+1;
    
end

ret(x,y)=tree(j).char-1;

end

%--------------------------------------------------------------------------
%Encoding Function
%-------------------------------------------------------------------------
function ret=encode(array,table)

ret=table(array(1,1)+1).b;
s=size(array);

for i=1:s(1)
    for j=1:s(2)
        if(~(i==1 && j==1))
            ret=[ret table(array(i,j)+1).b];
        end
    end
end

end

%--------------------------------------------------------------------------
%Traversing Function
%--------------------------------------------------------------------------
function s=traverse(tree)
s =struct('b',-1);

for i=1:256
    s(i).b=-1;
end

indx=length(tree);

for i=1:indx
    if tree(i).c1==-1 && tree(i).c2==-1
        j=i;
        steps=0;
        while j<=indx && tree(j).prnt~=-1
            if tree(tree(j).prnt).c1==j
                steps=[1 steps];
            else
                steps=[0 steps];
            end
            j= tree(j).prnt;
        end
        
        steps=steps(1:length(steps)-1);
        s(tree(i).char).b=steps;
        
    end
end


end
%--------------------------------------------------------------------------
%Huffman Tree Function
%--------------------------------------------------------------------------
%function output = function_name(input);
function tree=build_huffman(hist)

top=1;
indx=1;
j=1;
array =struct('c1', -1, 'c2',-1,'val',0,'char',-1,'i',-1,'prnt',-1);
tree =struct('c1', -1, 'c2',-1,'val',0,'char',-1,'i',-1,'prnt',-1);
x=0;
for i=1:256
    if hist(i)>0
        x=x+1;
        array(j).val=hist(i);
        array(j).char=i;
        array(j).c1=-1;
        array(j).c2=-1;
        array(j).i=-1;
        array(j).prnt=-1;
        j=j+1;
    end
end


while length(array) > 1
    
    [blah, order] = sort([array(:).val],'ascend');
    array = array(order);
    
    %check child 1
    if array(1).i==-1
        array(1).i=indx;
        tree(indx)=array(1);
        indx=indx+1;
    end
    
    %check child 2
    if array(2).i==-1
        array(2).i=indx;
        tree(indx)=array(2);
        indx=indx+1;
    end
    
    array(2).val=array(2).val+array(1).val;
    array(2).c1=array(1).i;
    array(2).c2=array(2).i;
    array(2).i=indx;
    array(2).char=-1;
    tree(array(2).c1).prnt=indx;
    tree(array(2).c2).prnt=indx;
    
    tree(indx)=array(2);
    indx=indx+1;
    
    array=array(2:length(array));
end

indx=indx-1;
tree(indx).prnt=-1;

end
%--------------------------------------------------------------------------