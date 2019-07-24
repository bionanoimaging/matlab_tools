% [myerr,myder] = matMulVec(aMatrix, aVector) : Multiplies a Vector (anImg) with a Matrix (Multfoectors). AnImg is always interpreted as a vector. The outermost dimension of MulFactors forms the result axis direction.
% In Einstein notation: result_k = aMatrix_nk * aVector_n
% It also works for the following case: result_k = aMatrix_nmk * aVector_nm
function res=matMulVec(aMatrix, aVector)
resDip=0
if isDipImage(aMatrix)
    aMatrix = double(aMatrix)
    resDip=1
end
if isDipImage(aVector)
    aVector = double(aVector)
    resDip=1
end
sm = size(aMatrix)
if ndims(aMatrix) > 2
    aMatrix = reshape(aMatrix,[prod(sm(1:end-1)) sm(end)])
end
sv = size(aVector)

aVector = reshape(aVector,[prod(sv),1])

res = aMatrix * aVector
if resDip
    res = dip_image(res)
end
