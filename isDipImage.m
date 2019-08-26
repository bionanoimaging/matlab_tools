% res = isDipImage(anImg): checks whether the input is of type DipImage independent whether cudaMat is used or not
function res = isDipImage(anImg)
if isa(anImg,'cuda')
    res=isDip(anImg);
else
    res=isa(anImg,'dip_image');
end
