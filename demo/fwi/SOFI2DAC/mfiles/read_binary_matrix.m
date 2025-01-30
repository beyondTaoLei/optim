function matrix=read_binary_matrix(nx,ny,matrix_file);

[fid,mess]=fopen(matrix_file,'r');
matrix_dummy=fread(fid,'float32');
matrix=reshape(matrix_dummy,ny,nx);
fclose(fid);

end
