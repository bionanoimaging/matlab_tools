function data=DimEdges(data)
global para;
for ind=1:length(data)
  if para.dim_edgeMethod ==2
    data{ind}=kwindow(data{ind},{'window',4;'mirror',1;'relXi',para.dim_above;'relXo',1.0+(1-para.dim_above);'relYi',para.dim_above;'relYo',1.0+(1-para.dim_above);'relZi',1.1;'relZo',1.2;'rect',1});
  elseif para.dim_edgeMethod ==3
      [dat,noiseRed]=DampEdge(data{ind},para.dim_above,3,1,2,para.dim_sigma);
      data{ind}=dat;
      para.noiseRed{ind}=noiseRed;
     % error('dim_edgeMethod 3 not implemented yet.');
  elseif para.dim_edgeMethod ==1
    data{ind}=kwindow(data{ind},{'window',4;'relXi',para.dim_above;'relXo',1.0;'relYi',para.dim_above;'relYo',1.0;'relZi',1.1;'relZo',1.2;'rect',para.dim_edgeShape});
  else
      error('Unknown dim_edgeMethod value. Choose 1 to 3');
  end
end
