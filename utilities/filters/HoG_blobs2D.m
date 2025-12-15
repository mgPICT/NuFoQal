function [outIm,whatScale,Direction] = HoG_blobs2D(I, options)
% This function HoG_blobs2D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to blobs, based on
% the method described by Frangi:2001 (Chapter 2) and its code.
%
% [J,Scale,Direction] = HoG_blobs2D(I, Options)
%
% inputs,
%   I : The input image (vessel image)
%   Options : Struct with input options,
%       .HoGScaleRange : The range of sigmas used, default [1 10]
%       .HoGScaleRatio : Step size between sigmas, default 2
%       .HoGAlpha : Frangi correction constant, default 0.5
%       .HoGBeta : Frangi correction constant, default 0.5
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)   
%
% Example,
%   I=double(imread ('vessel.png'));
%   Ivessel=FrangiFilter2D(I);
%   figure,
%   subplot(1,2,1), imshow(I,[]);
%   subplot(1,2,2), imshow(Ivessel,[0 0.25]);
%
% Written by Marc Schrijver, 2/11/2001
% Re-Written by D.Kroon University of Twente (May 2009)

defaultoptions = struct('HoGScaleRange', [1 10], 'HoGScaleRatio', 2, 'HoGAlpha', 0.5, 'HoGBeta', 0.5, 'verbose',true,'BlackWhite',true);

% Process inputs
if(~exist('options','var'))
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options)))
        warning('HoG_blobs2D:unknownoption','unknown options found');
    end
end

sigmas=options.HoGScaleRange(1):options.HoGScaleRatio:options.HoGScaleRange(2);
sigmas = sort(sigmas, 'ascend');

% Make matrices to store all filterd images
ALLfiltered=zeros([size(I) length(sigmas)]);
ALLangles=zeros([size(I) length(sigmas)]);

% Frangi filter for all sigmas
for i = 1:length(sigmas)
    % Show progress
    if(options.verbose)
        disp(['Current HoG Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Make 2D hessian
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));
    
    % Correct for scale
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
   
    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda2,Lambda1,Ix,Iy]=eig2image(Dxx,Dxy,Dyy);

    % Compute the direction of the minor eigenvector
    angles = atan2(Ix,Iy);
   
    % Compute the output image
    % The blobness Features
    if options.BlackWhite
        tra = Lambda1+Lambda2;
        det = imcomplement(Lambda1.*Lambda2);
    else
        tra = imcomplement(Lambda1+Lambda2);
        det = Lambda1.*Lambda2;
    end
    tra(tra<0) = 0;
    det(det<0) = 0;
	
    % Free memory
    clear LambdaAbs1 LambdaAbs2

    %Compute blobness function
    tra = tra./2;
    det = det.^(1/2);
    
    A = range(tra(:)).*options.HoGAlpha;
    B = range(det(:)).*options.HoGBeta;
    
    expTra = (1-exp(-(tra./A)));
    expDet = (1-exp(-(det./B)));
    
    Ifiltered = expTra.* expDet;
    
    % see pp. 45
    if(options.BlackWhite)
        Ifiltered(Lambda1<0)=0;
    else
        Ifiltered(Lambda1>0)=0;
    end
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
end

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1
    [outIm,whatScale] = max(ALLfiltered,[],3);
    outIm = reshape(outIm,size(I));
    if(nargout>1)
        whatScale = reshape(whatScale,size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
    end
else
    outIm = reshape(ALLfiltered,size(I));
    if(nargout>1)
            whatScale = ones(size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles,size(I));
    end
end
