function addXYtofeatures( input_feat_file, input_image, out_file )
%ADDXYTOFEATURES Takes a farsight file with the image features and a
%labeled image and adds the XY coordinates to the file
fid   = fopen( input_feat_file,'r' );
heads = fgetl( fid );
%tempA = fscanf( fid, '%g' );
A  = dlmread(input_feat_file,'\t');% reshape(tempA, [], numf );
numf  = size(A,2);
fprintf( 1, 'Assuming the number of features are %d \n',numf);

im = imread( input_image );
max_lab = max(max(im));
ind = 1;
centroid = zeros(size(A,1),2);
for i=1:max_lab
	[x,y] = find(im==i);
	if ~isempty(x)
		centroid(ind,1) = round(mean(x));
		centroid(ind,2) = round(mean(y));
		ind=ind+1;
	end
	x = [];
	y = [];
end
new_feat = zeros(size(A,1),4);
new_feat(:,1:2) = centroid;
new_feat(:,3)   = A(:,2);
new_feat(:,4)   = A(:,numf-1);
 A(1,:)
%fod = fopen( out_file, 'w' );
%fwrite(fod, heads );
%fprintf(fod,'%.2f\t%.2f\t%.2f\t%.2f\n',new_feat);
dlmwrite(out_file,new_feat,'\t');
end
