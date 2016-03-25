function width0=PIPE_estimate_width0(data,coords)
    %%% This function receives a 3D data matrix and the center-of-pulse
    %%% coordinates and estimates the width of the initial photo-converted
    %%% intensity profile.
    
    %%% unpack coordinates
    temp=num2cell(coords);
    [t0,y0,x0,direction,thickness]=temp{:};
    
    %%% define and get initial intensity profile
    [profile_mask,perp_direction]=PIPE_define_1d_profile_mask(y0,x0,direction,thickness,size(data(:,:,1)));
    profile=PIPE_apply_profile_mask(data(:,:,t0),profile_mask, perp_direction, thickness);
    if isempty(profile)
        width0=0;
        return
    end
    
    %%% fit profile to Gaussian
    [a0,b0]=max(profile);
    c0=numel(profile)/4;
    d0=0;
    gauss=@(a,b,c,d,x)a*exp(-((x-b)/c).^2)+d;
    [model,params]=fit((1:numel(profile))',profile(1:numel(profile)),gauss,'Startpoint',[a0 b0 c0 d0],'Lower',[0 1 0 -max(profile)],'Upper',[100*(max(profile)-min(profile)) 100*numel(profile) 100*numel(profile) max(profile)]);
    if params.rsquare>0.65
        width0=model.c;
    else
        width0=0;
    end
end