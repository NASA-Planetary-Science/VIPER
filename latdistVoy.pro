FUNCTION latdistVoy, scpos, pos, no, tpar, mass, char, fcorot, LATDIST=LATDIST
  
  ; Inputs
  ; 
  ; 
  ; 
  ; Keyword Arguments
  ; 
  ; Outputs
  
  ; Constants
  degrad=!PI/180.
  te = 50d ; May need to change this Logan 10/17/16
  
  ; Create arrays of correct size
  ns     = N_ELEMENTS(no)
  ni     = FLTARR(ns+1)
  ni(ns) = TOTAL(no[*]*char[*]) ; Total electron density due to quasineutrality
  FOR i = 0, ns-1 DO ni(i) = no(i)*ni(ns)

  ; Do a distribution along the field line if requested
  IF KEYWORD_SET(LATDIST) THEN BEGIN
    ; If this keyword is set, then the scpos is assumed to be where the density is measured, which if this function was
    ; used beforehand to get the density at the location of the centrifugal equator, then the location of the centrifugal
    ; equator will be the scpos of that measurement
    
    ;nl=#latitude bins, dl=interval(degrees)
    nl = 85
    dl = 1.0
    
    pos_r   = FLTARR(nl)
    pos_l   = FLTARR(nl)
    lat     = FLTARR(nl)
    r       = FLTARR(nl)
    ngrid   = FLTARR(nl,ns+1)
    ngrid2  = FLTARR(nl,ns+3)
    ngrid(*) = 0.
    
    ; pos(0) = r   - radial distance
    ; pos(1) = lat - centrifugal latitude
    pos = FLTARR(2)
    
    ; FIND POSITION WHERE DENSITY TO BE DETERMINED
    FOR l = 0, nl-1 DO BEGIN
      lat(l) = l*dl
      pos(1) = lat(l) + scpos(1) ; Augment by actual location so that it steps away from measurement location LOGAN 10/17/16
      ; Assumption: Dipole and magnetic and rotational equators aligned
      fact   = cos(pos(1)*degrad)/cos(scpos(1)*degrad)
      pos(0) = scpos(0)*fact*fact
      pos_r(l) = pos(0)
      pos_l(l) = pos(1)
      
      FORWARD_FUNCTION BBeq
      BB = BBeq(pos(0),lat(1))
      BeB = alog(BB)
      
      anis    = MAKE_ARRAY(N_ELEMENTS(ni))
      anis[*] = 1d
      anis[-1] = 2d
      
      no = [no[*], ni[-1]]
      ; MULTIPLE MAXWELLIANS - isotropic
      FORWARD_FUNCTION nmax
      ngrid(l,*) = nmax(ns,pos,scpos,no,tpar,anis,BeB,te,mass,char,fcorot)
      ;ngrid(l,*) = nmax_Voy(ns,pos,scpos,ni,tpar,te,mass,char,fcorot)
      
      
    ENDFOR
    
    ; Fill array that includes latitude and position, densities start on 3rd instead of first
    FOR ii = 0, ns DO BEGIN
      ngrid2(*,ii+2) = ngrid(*,ii)
    ENDFOR
    ngrid2(*,0)      = pos_r(*)
    ngrid2(*,1)      = pos_l(*)
    
  ENDIF ELSE BEGIN
    
    FORWARD_FUNCTION BBeq   
    BB = BBeq(pos(0),pos(1))
    BeB = alog(BB)

    anis    = MAKE_ARRAY(N_ELEMENTS(ni))
    anis[*] = 1d
    anis[-1] = 2d
    
    FORWARD_FUNCTION nmax
    no = [no[*], ni[-1]] ; Dougherty 2020_08_21 - I think we were squaring the density basically before
    ngrid = nmax(ns,pos,scpos,no,tpar,anis,BeB,te,mass,char,fcorot)
    ;ngrid = nmax_Voy(ns,pos,scpos,ni,tpar,te,mass,char,fcorot)   
 
  ENDELSE
  IF KEYWORD_SET(LATDIST) THEN RETURN, ngrid2 ELSE RETURN, ngrid
END
