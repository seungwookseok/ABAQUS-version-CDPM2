c Concrete damage-plasticity model 2 (CDPM2) _ VUMAT for ABAQUS
c
c The CDPM2 was orgiginally developed by the research group of Dr. Peter Grassl (Univ. of Glasgow, UK)
c and has been implemented in LS-DYNA as MAT CDPM (MAT 273).
c Afterwards, it was revised for use in ABAQUS by Seungwook Seok (a PhD student in Civil Eng. at Purdue Univ., USA)
c updated in 2019-07-21
c
c*****************************************************************
c Support page: http://petergrassl.com/Research/DamagePlasticity/CDPMLSDYNA/index.html
c               https://github.com/seungwookseok/ABAQUS-version-CDPM2
c
c Key references: 
c 1) P. Grassl, D. Xenos, U. Nyström, R. Rempling, K. Gylltoft.: "CDPM2: A damage-plasticity approach to modelling the failure of concrete". International Journal of Solids and Structures. Volume 50, Issue 24, pp. 3805-3816, 2013.
c 2) P. Grassl, U. Nyström, R. Rempling and K. Gylltoft, "A damage-plasticity model for the dynamic failure of concrete", 8th International Conference on Structural Dynamics, Leuven, Belgium, 2011
c 3) P. Grassl and M. Jirasek: "Damage-plastic model for concrete failure". International Journal of Solids and Structures. Vol. 43, pp. 7166-7196, 2006.
c*****************************************************************
c
c-----------------------------------------------------------------
c # of properties (props) = 13
c # of state variables (state) = 49
c-----------------------------------------------------------------
c User-defined material properties are as follows:
c
c props(1) unitflag: Units flag (US Costumary units = 0, SI units = 1)
c props(2) ym: Young’s modulus
c props(3) pr: Poisson’s ratio
c props(4) fc': Concrete uniaxial compressive strength
c props(5) ft: Concrete uniaxial tensile strength
c props(6) gft: Fracture energy 
c props(7) gfc: Crushing energy (not used, to be updated)
c props(8) ah: Parameter of the hardening ductility measure
c props(9) bh: Parameter of the hardening ductility measure
c props(10) ch: Parameter of the hardening ductility measure
c props(11) dh: Parameter of the hardening ductility measure
c props(12) hp: Hardening modulus for qh2
c props(13) as: Parameter of the softening ductility measure
c-----------------------------------------------------------------
c The state variables are stored as:
c
c state(*,1) ---------- kappa
c state(*,2) ---------- equivalent strain
c state(*,3) ---------- plastic strain xx (or 11)
c state(*,4) ---------- plastic strain yy (or 22)
c state(*,5) ---------- plastic strain zz (or 33)
c state(*,6) ---------- plastic strain xy (or 12)
c state(*,7) ---------- plastic strain yz (or 23)
c state(*,8) ---------- plastic strain xz (or 13)
c state(*,9) ---------- kappa tension kdt
c state(*,10) --------- kappa tension 1 kdt1
c state(*,11) --------- kappa tension 2 kdt2
c state(*,12) --------- kappa compression kdc
c state(*,13) --------- kappa compression 1 kdc1
c state(*,14) --------- kappa compression 2 kdc2
c state(*,15) --------- damage variable tension omegaT (wt)
c state(*,16) --------- damage variable tension omegaC (wc)
c state(*,17) --------- strain rate factor used for the dynamic formulation based on quasi-static analysis
c state(*,18) --------- alphac is the compression factor given in paper in IJSS by Grassl et al. in equation 46
c state(*,19) --------- equivalent strain tension eqstrT
c state(*,20) --------- equivalent strain compression eqstrC
c state(*,21) --------- total strain along xx (or 11)
c state(*,22) --------- total strain along yy (or 22)
c state(*,23) --------- total strain along zz (or 33)
c state(*,24) --------- total strain along xy (or 12)
c state(*,25) --------- total strain along yz (or 23)
c state(*,26) --------- total strain along xz (or 13)
c state(*,27) --------- equivalent strain (without rate factor influence)
c state(*,28) --------- element deletion flag: "1" active, "0" inactive (deleted element)
c state(*,29:49) ------ to be updated
c-----------------------------------------------------------------

      subroutine vumat(
c Read only (unmodifiable)variables -
     &  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     &  stepTime, totalTime, dt, cmname, coordMp, charLength,
     &  props, density, strainInc, relSpinInc,
     &  tempOld, stretchOld, defgradOld, fieldOld,
     &  stressOld, stateOld, enerInternOld, enerInelasOld,
     &  tempNew, stretchNew, defgradNew, fieldNew,
c Write only (modifiable) variables -
     &  stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     &  charLength(nblock), strainInc(nblock,ndir+nshr),
     &  relSpinInc(nblock,nshr), tempOld(nblock),
     &  stretchOld(nblock,ndir+nshr),
     &  defgradOld(nblock,ndir+nshr+nshr),
     &  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     &  stateOld(nblock,nstatev), enerInternOld(nblock),
     &  enerInelasOld(nblock), tempNew(nblock),
     &  stretchNew(nblock,ndir+nshr),
     &  defgradNew(nblock,ndir+nshr+nshr),
     &  fieldNew(nblock,nfieldv),
     &  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     &  enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname

      real*8 eps(6);
      integer maxnip,lft,llt
c     maxnip    ---- variable denoting max number of integration points
      integer mx,i,j

      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

c     ym   --------------- Young's modulus
c     pr   --------------- Poisson's ratio
c     ecc  --------------- eccentricity parameter used in the plasticity law. Its calibration is described in Jirasek & Bazant (2002)
c     qh0  --------------- value of the derivative of the 1st hardening function at kappa=0 as described in eq. (30) of the IJSS paper by Grassl et al.
c     ft   --------------- tensile strength
c     fc   --------------- compressive strength
c     hp   --------------- value of the derivative of the 2nd hardening function as described in eq. (31) of the IJSS paper by Grassl et al.
c     ah,bh,ch,dh -------- hardening parameters used in eq. (33) of the IJSS paper by Grassl et al.
c     as,bs,df ----------- softening  parameters used in eqs. (56) and (50) of the IJSS paper by Grassl et al. 
c     fc0  --------------- reference compressive strength used in the definition of the rate factor
c     type --------------- softening type used in the formulation of the tensile damage variable
c     wf   --------------- max crack opening displacement used in the tensile damage law
c     wf1  --------------- max crack opening displacement used in the bilinear tensile damage law
c     efc  --------------- compressive strain used as a parameter in the exponential compressive damage law
c     ft1  --------------- tensile stress used in the billinear tensile damage law
c     strrateflg --------- variable denoting whether impact phenomena are taken into account in the constitutive law
c     failflg ------------ flag denoting when an element should be deleted
c     m0   --------------- parameter used in the plasticity law calculated in eq.(20) of the IJSS paper by Grassl et al.
c     isoflag------------- flag to denote whether isotropic or anisotropic damage law is used

c Variables used for the evaluation of the plasticity algorithm
      real*8 sig(6),sigEff(6),oldStrain(6),convStrain(6),
     $     deltaTotStrain(6),
     $     tempTotStrain(6),elStrain(6),princStress(3),
     $     princDir(3,3),totStrain(6),plastStrain(6),
     $     sum,tempTheta,
     $     strain(6), sigVTrial,rhoTrial,thetaTrial,apexStress,
     $     tempkappaP,yieldval
      integer subincCounter,rtype,subincFlag,converged,l,k
c       totstrain      -----------  total strain vector equal to the sum of  plastic and elastic strains
c       sigVTrial      ----------- trial volumetric stress 
c       rhoTrial       ----------- trial deviatoric stress
c       thetaTrial     ----------- trial Lode angle
c       apexStress     ----------- apexstress. Used only in the vertex case of the plasticity algorithm
c       yieldval   - value of the yield function
c       subincCounter   - counter of subincrementations performed in the plasticity algorithm
c       subincFlag   - flag denoting whether or not subincrementation is taking place
c                      =0 no subincrementation is taking place
c                      =1 subincrementation is taking place
c       rtype      - return type of the yield surface
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged       - integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     sigEff(6)        ---------- effective stress (no damage included)
c     oldStrain(6)     ---------- old elastic strain 
c     convStrain(6)    ---------- last elastic strain vector for which the plasticity algorithm has converged
c     deltaTotStrain(6)   ---------- difference between last and new strain vector
c     tempTotStrain(6)  ---------- temporary elastic strain
c     elStrain(6)      ----------  elastic strain
c     princStress(3)   ----------  array containing principal stresses
c     princDir(3,3)    ----------  matrix containing principal eigenvectors stored columnwise
c     plastStrain(6)   ----------  plastic strain vector
c     sum              ----------  variable used for summation
c     tempTheta        ---------- temporary Lode angle
c     strain(6)        ---------- strain array
c     l,k              ---------- counters used in various iterations
c
c Variables used in the damage algorithm
      real*8 alpha,omegaT,omegaC,effStressT(6),effStressC(6),
     $     rateFactor,cdpm2u_computeRateFactor,strainrate(6),epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     length,deltaElStrain(6),sigOld(6), damagedGP    
c     alpha             --------- variable used to identify amount of contribution of compressive stresses and is defined in eq. (46) of the IJSS paper by Grassl et al.
c     omegaT            --------- tensile damage variable
c     omegaC            --------- compressive damage variable
c     effStressT        --------- effective tensile part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     effStressC        --------- effective compressive part of the effective stress tensor(no damage included) as described in Section 2.1 of the IJSS paper by Grassl et al.
c     rateFactor       ---------- variable used to incorporate impact effects on the constitutive law
c     cdpm2u_computeRateFactor --------- function used to calculate the rate factor
c     strainrate(6)    ---------- rate of the strain tensor used for the calculation of the rateFacto
c     epsilonT         ---------- tensile equivalent strain
c     epsilonC         ---------- compressive equivalent strain
c     kappaDT          ---------- history parameter kappaDT        
c     kappaDT1         ---------- history parameter kappaDT1       
c     kappaDT2         ---------- history parameter kappaDT2       
c     kappaDC          ---------- history parameter kappaDC        
c     kappaDC1         ---------- history parameter kappaDC1       
c     kappaDC2         ---------- history parameter kappaDC2       
c     length           ---------- characteristic length used for the formulation of the tensile damage function based on the crack-band approach
c     deltaElStrain(6) ---------- elastic strains of the previous step
c     damagedGP        ---------- number of damaged gausspoints within analysed element 
c
c-----------------------------------------------------------------
c
c Input variables
      real*8 failflag
c     type (Damage type)
c		     = 0.0: Linear softening     
c		     = 1.0: Bi-Linear softening  
c		     = 2.0: Exponential softening
c		     = 3.0: NO DAMAGE
c     strrateflg (Strain rate parameter)
c		     = 0.0: No rate effects
c		     = 1.0: Rate effects included  
c     failflag 
c		     = 0.0: if not ALL  gausspoints of the element are damaged        
c            = 1.0: if ALL gausspoints of the element are damaged  
c     unitflag
c            = 0.0: US Costumary units
c            = 1.0: SI units
c
c Additional variables
      real*8 fbfcRatio,fb,epsforecc,maxsubinc
      real*8 jacobian(4,4),residuals(4),normResiduals,jacobianPrev(4,4),
     $       jacobianOld(4,4),residualsOld(4),normResidualsOld,
     $       endflag,jacDet,jacDetOld,
     $       evecOld(4,4),evalOld(4),itNoOld,roNoOld,
     $       evec(4,4),eval(4),itNo,roNo
c
c-----------------------------------------------------------------

c Default parameters:

      type = 1.0
      strrateflg = 0.0
      failflg = 1 ! 0.0: not active, x>0.0: active and element will erode to 1
                    ! in x percent of the integration points. If x=0.60, 60% of
                    ! all integration points must fail before erosion.
      maxnip = 1 ! number of the integration points for C3D8R
      isoflag = 0.0
      bs = 1.0
      df = 0.85
      printflag = 0.0 ! means nothing!
      maxsubinc = 20
      endflag = 0

c Material constants

      unitflag = props(1)
      ym = props(2)
      pr = props(3)
      fc = props(4)
      ft = props(5)
      gft = props(6)
c      gfc = props(7)
c      gfc = 25 ! (N/mm)
      ah = props(8)
      bh = props(9)
c      bh = 0.0025 ! Bh
c      bh = -2.29*pepeak + 0.00046 ! (A.10), Grassl and Jirasek, p.25
      ch = props(10)
      dh = props(11)
      hp = props(12)
      as = props(13)
	  
      
c This is the global tolerance. All other tolerances are made relative to this global tolerance

      gTol = 1.e-4
      damagedGP = 0.
	  
      if (ym<0 ) then
         isoflag = 1.0
         ym = abs(ym)
      end if
	  	  
c Calculated parameters
      if (unitflag == 0.0) then
         unitconv = 6.89475908677537 ! Convert ksi to MPa
      else if (unitflag == 1.0) then
         unitconv = 1
      end If
      fbfcRatio = 1.5*(fc*unitconv)**(-0.075) ! (10), Papanikolaou and Kappos (2007), p.5
      fb = fbfcRatio*fc
c      fb = 1.16*fc ! Kupfer et al. (1969)
      fc0 = ((fc*unitconv)**(1.855)/60.) / unitconv ! (14), Papanikolaou and Kappos (2007), p.8	  
	  qh0 = fc0/fc
      epsforecc = ft/fb * (fb**2.-fc**2.)/(fc**2.-ft**2.)
      ecc = (1.+epsforecc) / (2.-epsforecc)
      m0 = 3. * (fc**2.-ft**2.)/(fc*ft) * ecc/(ecc+1.)
      e0 = ft/ym ! e0
      wf = Gft/(0.225*ft) ! see p.7 in the paper
      wf1 = 0.15*wf ! Jirasek and Zimmermann (1998)
      ft1 = 0.3*ft ! Jirasek and Zimmermann (1998)

      lft = 1 ! number of first material point?
      llt = nblock ! number of last material point?
	  
      if (stepTime .eq. 0. .and. lanneal .eq. 0) then
         call cdpm2u_printallinputvariables()
      end if
	  
      do i = lft,llt
	  
         stateNew(i,28) = stateOld(i,28)
         tempkappaP = stateOld(i,1)
         length = charLength(i)
		 
         efc = 0.0001
	 
         jacobianOld(1,1) = stateOld(i,29)
         jacobianOld(1,2) = stateOld(i,30)
         jacobianOld(1,3) = stateOld(i,31)
         jacobianOld(1,4) = stateOld(i,32)
         jacobianOld(2,1) = stateOld(i,33)
         jacobianOld(2,2) = stateOld(i,34)
         jacobianOld(2,3) = stateOld(i,35)
         jacobianOld(2,4) = stateOld(i,36)
         jacobianOld(3,1) = stateOld(i,37)
         jacobianOld(3,2) = stateOld(i,38)
         jacobianOld(3,3) = stateOld(i,39)
         jacobianOld(3,4) = stateOld(i,40)
         jacobianOld(4,1) = stateOld(i,41)
         jacobianOld(4,2) = stateOld(i,42)
         jacobianOld(4,3) = stateOld(i,43)
         jacobianOld(4,4) = stateOld(i,44)
         residualsOld(1) = stateOld(i,45)
         residualsOld(2) = stateOld(i,46)
         residualsOld(3) = stateOld(i,47)
         residualsOld(4) = stateOld(i,48)
         normResidualsOld = stateOld(i,49)
		 
c Write strain rate vector

         eps(1:6) = strainInc(i,1:6)

         do l = 1,6
            totStrain(l) = eps(l) + stateOld(i,l+20)
            strainrate(l) = eps(l)
            plastStrain(l) = stateOld(i,l+2)
            oldStrain(l) = stateOld(i,l+20)
            convStrain(l) = oldStrain(l)
            tempTotStrain(l) = totStrain(l)
            deltaTotStrain(l) = eps(l)
         enddo

         subincCounter=0
         subincFlag=0
         
         converged=1
                 
         do while ( converged .eq. 1 .or. subincFlag .eq. 1 ) 
            
            do l=1,6
               elStrain(l) = tempTotstrain(l) - plastStrain(l)              
            enddo   

            call cdpm2u_computeStressesfromStrains(sigEff,elStrain,
     $           ym,pr)
            
            call cdpm2u_computeTrialCoordinates(sigEff,sigVTrial,
     $           rhoTrial,tempTheta)
            thetaTrial=tempTheta
            call cdpm2u_computeYieldValue(yieldval,sigVTrial,rhoTrial,
     $           thetaTrial,tempKappaP)

            apexStress=0.

            if (yieldval .gt. 0.) then
               call cdpm2u_checkForVertexCase(apexStress,sigVTrial,
     $              tempKappaP,rtype)
               if (rtype.eq.1 .or. rtype .eq. 2) then
                  call cdpm2u_performVertexReturn(sigEff,
     $                 apexStress,tempKappaP,rtype,converged)
               end if
               
               if (rtype.eq.0) then
                  call cdpm2u_performRegularReturn(sigEff, 
     $                 tempKappaP,converged,ym,pr,gTol,jacobian,
     $                 residuals,normResiduals,endflag)
               end if
            else
               converged = 0
               do l=1,6        
                  plastStrain(l) = stateOld(i,l+2)
               enddo
               goto 925
            end if                     
            
            if ( converged .eq. 1 ) then
			   
               subincCounter=subincCounter+1
               if ( subincCounter .gt. maxsubinc ) then
                  call cdpm2u_computeMatrixDeterminant(4,jacobian,jacDet)
                  call cdpm2u_computeMatrixDeterminant(4,jacobianOld,jacDetOld)
                  write(*,*) '*** Perform Plasticity return with' 
                  write(*,*) 'subincrementation methodology'
                  write(*,*) 'No convergence reached !***'
                  write(*,*) 'jacobianOld(1,:)', jacobianOld(1,:)
                  write(*,*) 'jacobianOld(2,:)', jacobianOld(2,:)
                  write(*,*) 'jacobianOld(3,:)', jacobianOld(3,:)
                  write(*,*) 'jacobianOld(4,:)', jacobianOld(4,:)
                  write(*,*) 'residualsOld', residualsOld
                  write(*,*) 'normResidualsOld', normResidualsOld
                  write(*,*) 'jacobianPrev(1,:)', jacobianPrev(1,:)
                  write(*,*) 'jacobianPrev(2,:)', jacobianPrev(2,:)
                  write(*,*) 'jacobianPrev(3,:)', jacobianPrev(3,:)
                  write(*,*) 'jacobianPrev(4,:)', jacobianPrev(4,:)
                  write(*,*) 'jacobian(1,:)', jacobian(1,:)
                  write(*,*) 'jacobian(2,:)', jacobian(2,:)
                  write(*,*) 'jacobian(3,:)', jacobian(3,:)
                  write(*,*) 'jacobian(4,:)', jacobian(4,:)
                  write(*,*) 'jacDet', jacDet
                  write(*,*) 'jacDetOld', jacDetOld
                  write(*,*) 'residuals', residuals
                  write(*,*) 'normResiduals', normResiduals
                  call jacobi_eigenvalue (4,jacobianOld,1000,evecOld,
     $                 evalOld,itNoOld,roNoOld)
                  call jacobi_eigenvalue (4,jacobian,1000,evec,
     $                 eval,itNo,roNo)
                  write(*,*) 'evalOld', evalOld
                  write(*,*) 'eval', eval
                  endflag = 1
                  call cdpm2u_performRegularReturn(sigEff, 
     $                 tempKappaP,converged,ym,pr,gTol,jacobian,
     $                 residuals,normResiduals,endflag)
                  stop
               else if (subincCounter .gt. maxsubinc-1 .and. 
     $                 tempKappaP .lt. 1.0 ) then
                  tempKappaP = 1.
               end if
               subIncFlag = 1
               do l = 1,6
                  deltaTotStrain(l) = deltaTotStrain(l)*0.5
                  tempTotStrain(l) = convStrain(l)+deltaTotStrain(l)
               enddo
			   
               jacobianPrev = jacobian
			   
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 0) then
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $          ym,pr)
               do l=1,6        
                  plastStrain(l) = totStrain(l)-elStrain(l)
               enddo
            else if ( converged .eq. 0 .and. 
     $              subIncFlag .eq. 1) then
c               write(*,*) '*** Subincrementation required',subincCounter               
               call cdpm2u_computeStrainsfromStresses(sigEff,elStrain,
     $              ym,pr)
               do l = 1,6
                  plastStrain(l) = tempTotStrain(l)-elStrain(l)
                  convStrain(l) = tempTotStrain(l)
                  deltaTotStrain(l) = totStrain(l)-convStrain(l)
                  tempTotStrain(l) = totStrain(l)
               enddo
               subincCounter = 0
               subincFlag = 0
               converged = 1            
            end if
         end do
         
                     
 925     continue
         if (type .eq. 3.0) then
            omegaT = 0.0
            omegaC = 0.0            
            epsilonT = 0.0
            kappaDT = 0.0
            kappaDT1 = 0.0
            kappaDT2 = 0.0
            kappaDC = 0.0
            kappaDC1 = 0.0
            kappaDC2 = 0.0
            omegaT = 0.0
            omegaC = 0.0
            rateFactor = 0.0
            alpha = 0.0
            epsilonT = 0.0
            epsilonC = 0.0
            do l = 1,6
               sig(l) = sigEff(l)
            enddo            
            goto 152
         end if

c Initialize parameters used in the damage algorithm             

         rateFactor = stateOld(i,17)
         epsilonT = stateOld(i,19)
         epsilonC = stateOld(i,20)
         kappaDT = stateOld(i,9)
         kappaDT1 = stateOld(i,10)
         kappaDT2 = stateOld(i,11)
         kappaDC = stateOld(i,12)
         kappaDC1 = stateOld(i,13)
         kappaDC2 = stateOld(i,14)
         omegaT = stateOld(i,15)
         omegaC = stateOld(i,16)
         alpha = stateOld(i,18)

         call cdpm2u_computeAlpha(effStressT,effStressC,sigEff,alpha)
         sum = 0.
         do l = 1,6
c        Compute norm of increment of plastic strains       
            sum = sum + (plastStrain(l)-stateOld(i,l+2))**2.
            deltaElStrain(l) = (stateOld(i,l+20)-stateOld(i,l+2))
         enddo
         call cdpm2u_computeStressesfromStrains(sigOld,deltaElStrain,
     $        ym,pr)
         sum = sqrt(sum)
         call cdpm2u_computeDamage(omegaC,omegaT,strainrate,
     $        rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,
     $        kappaDT2,kappaDC,kappaDC1,kappaDC2,sigEff,sum,
     $        tempKappaP,length,sigOld,
     $        stateOld(i,18),stateOld(i,27),stateNew(i,27))
		 
         do l = 1,6
            if (isoflag .eq. 1.0) then
               sig(l) = (1.-omegaT)*sigEff(l)
            else
               sig(l) = (1.-omegaT)*effStressT(l)
     $                  +(1.-omegaC)*effStressC(l)
            end if
         end do

 152     continue
         
c Write the history variable at the end of the routine

         stateNew(i,1) = tempkappaP 
         stateNew(i,2) = epsilonT
         stateNew(i,3) = plastStrain(1)
         stateNew(i,4) = plastStrain(2)
         stateNew(i,5) = plastStrain(3)
         stateNew(i,6) = plastStrain(4)
         stateNew(i,7) = plastStrain(5)
         stateNew(i,8) = plastStrain(6)
         stateNew(i,9) = kappaDT
         stateNew(i,10) = kappaDT1
         stateNew(i,11) = kappaDT2
         stateNew(i,12) = kappaDC
         stateNew(i,13) = kappaDC1
         stateNew(i,14) = kappaDC2
         stateNew(i,15) = omegaT
         stateNew(i,16) = omegaC
         stateNew(i,17) = rateFactor
         stateNew(i,18) = alpha
         stateNew(i,19) = epsilonT
         stateNew(i,20) = epsilonC
         stateNew(i,21) = totStrain(1)
         stateNew(i,22) = totStrain(2)
         stateNew(i,23) = totStrain(3)
         stateNew(i,24) = totStrain(4)
         stateNew(i,25) = totStrain(5)
         stateNew(i,26) = totStrain(6)
		 
         stateNew(i,29) = jacobian(1,1)
         stateNew(i,30) = jacobian(1,2)
         stateNew(i,31) = jacobian(1,3)
         stateNew(i,32) = jacobian(1,4)
         stateNew(i,33) = jacobian(2,1)
         stateNew(i,34) = jacobian(2,2)
         stateNew(i,35) = jacobian(2,3)
         stateNew(i,36) = jacobian(2,4)
         stateNew(i,37) = jacobian(3,1)
         stateNew(i,38) = jacobian(3,2)
         stateNew(i,39) = jacobian(3,3)
         stateNew(i,40) = jacobian(3,4)
         stateNew(i,41) = jacobian(4,1)
         stateNew(i,42) = jacobian(4,2)
         stateNew(i,43) = jacobian(4,3)
         stateNew(i,44) = jacobian(4,4)
         stateNew(i,45) = residuals(1)
         stateNew(i,46) = residuals(2)
         stateNew(i,47) = residuals(3)
         stateNew(i,48) = residuals(4)
         stateNew(i,49) = normResiduals

		 
c Element deletion

c         if ( stateNew(i,15).gt. 0.9995 .and. stateNew(i,16).gt. 0.9995 ) then
c            damagedGP = damagedGP+1.0
c         end if
c         if (failflg .gt. 0.0) then
c            damagedGP = damagedGP/FLOAT(maxnip) ! e.g., maxnip = 1 for C3D8R
c            if (damagedGP .ge. failflg) then
c               stateNew(i,28) = 0
c            end if
c         end if
		 
         
c Write sig components

         stressNew(i,1:6) = sig(1:6)
		 
      end do
	  
      return
      end
	  
      
      subroutine cdpm2u_printallinputvariables()
c     Subroutine to check if all model parameters have values. If a parameter does not have any values a default value is provided.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag      
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

c Write output
c     Print all input variables
c
      write(*,*) '*****************************************************'
      write(*,*) '*****************************************************'
      write(*,*) '  ABAQUS is using the usermaterial CDPM2             '
      write(*,*) '  -------------------------------------------------  '
      write(*,*) '                                                     '
      write(*,*) '  Material parameters:                               '
      write(*,2) '   YM (Youngs modulus)............... = ',ym
      write(*,2) '   PR (Poissons ratio)............... = ',pr
      write(*,2) '   ECC (Eccentricity)................ = ',ecc
      write(*,2) '   QH0 (Initial hardening)........... = ',qh0
      write(*,2) '   FT (Uniaxial tension strength).... = ',ft
      write(*,2) '   FC (Uniaxial compression strength) = ',fc
      write(*,2) '   HP (Hardening modulus)............ = ',hp
      write(*,2) '   AH (Hardening ductility measure).. = ',ah
      write(*,2) '   BH (Hardening ductility measure).. = ',bh
      write(*,2) '   CH (Hardening ductility measure).. = ',ch
      write(*,2) '   DH (Hardening ductility measure).. = ',dh
      write(*,2) '   AS (Damage ductility measure)..... = ',as
      write(*,2) '   DF (Dilation constant)............ = ',df
      write(*,2) '   FC0 (Initial compression strength) = ',fc0
      write(*,2) '   TYPE (Damage type)................ = ',type
      write(*,*) '     EQ.0.0: Linear softening           '
      write(*,*) '     EQ.1.0: Bi-Linear softening        '
      write(*,*) '     EQ.2.0: Exponential softening      '
      write(*,*) '     EQ.2.0: No damage                  '
      write(*,2) '   BS (Damage: ductility parameter)   = ',bs
      write(*,2) '   WF (Damage: disp threshold 0)..... = ',wf
      write(*,2) '   WF1 (Damage: disp threshold 1).... = ',wf1
      write(*,2) '   FT1 (Damage: stress threshold 1).. = ',ft1
      write(*,2) '   EFC (strain threshold in comp).... = ',efc
      write(*,2) '   STRPR (Strain rate flag).......... = ',strrateflag
      write(*,*) '     EQ.0.0: No rate effects            '
      write(*,*) '     EQ.1.0: Rate effects included      '
      write(*,2) '   ISOFLAG (isotropic damage flag).. =  ',isoflag
      write(*,*) '     EQ.0.0: standard model with two damage params'
      write(*,*) '     EQ.1.0: model with one damage param  '

      write(*,*) '*****************************************************'
      write(*,*) '  History variables:                                 '
      write(*,*) '   Var  #1 = kappa'
      write(*,*) '   Var  #2 = equivalient strain'
      write(*,*) '   Var  #3 = plastic strain direction 1'
      write(*,*) '   Var  #4 = plastic strain dierction 2'
      write(*,*) '   Var  #5 = plastic strain direction 3'
      write(*,*) '   Var  #6 = hardening tension kdt'
      write(*,*) '   Var  #7 = hardening tension kdt1'
      write(*,*) '   Var  #8 = hardening tension kdt2'
      write(*,*) '   Var  #9 = hardening compression kdc'
      write(*,*) '   Var #10 = hardening compression kdc1'
      write(*,*) '   Var #11 = hardening compression kdc2'
      write(*,*) '   Var #12 = damage function tension wt'
      write(*,*) '   Var #13 = damage function compression wc'
      write(*,*) '   Var #14 = (internal book keeping)'
      write(*,*) '   Var #15 = compression factor alphac'
      write(*,*) '   Var #16 = ratefactor alphar'
      write(*,*) '   Var #17 = elastic strain direction 1'
      write(*,*) '   Var #18 = elastic strain direction 2'
      write(*,*) '   Var #19 = elastic strain direction 3'
      write(*,*) '   Var #20 = equivalent strain tension'
      write(*,*) '   Var #21 = equivalent strain compression'
      write(*,*) '   Var #22-#27 = undamaged stresses'
      write(*,*) '*****************************************************'
      write(*,*) '  Check that the bulk modulus and shear modulus are  '
      write(*,*) '  correctly calculated. They MUST be set on the      '
      write(*,*) '  input card.                                        '
      write(*,1) '  With E = ',ym,' and pr = ',pr
      write(*,2) '  BLK = E/(3*(1-2pr)) = ',ym/(3.0*(1.-2.*pr))
      write(*,2) '  SHR = E/(2*(1+pr))  = ',ym/(2.*(1.+pr))
      write(*,*) '*****************************************************'
      write(*,*) '*****************************************************'
 1        format(1x,A,1pE9.3,A,1pE9.3)
 2            format(1x,A,1pE12.5)

      return
      end    


      subroutine cdpm2u_computeStrainsfromStresses(stress,strain,ym,pr)
c     Subroutine to calculate elastic strains from the stress tensor. Performs operation epsilon = D : sigma
      real*8 stress(6),strain(6),ym,pr
c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      strain(1)=(stress(1) - pr * stress(2) - pr * stress(3))/ym
      strain(2)=(-pr*stress(1) + stress(2) - pr * stress(3))/ym
      strain(3)=(-pr*stress(1) - pr * stress(2) + stress(3))/ym
      strain(4)=(2. * (1+pr) * stress(4) )/ym
      strain(5)=(2. * (1+pr) * stress(5) )/ym
      strain(6)=(2. * (1+pr) * stress(6) )/ym
      return
      end
      
      subroutine cdpm2u_computeStressesfromStrains(stress,strain,ym,pr)
c     Subroutine to calculate strains from the elastic strain tensor. Performs operation sigma = C : epsilon

c     strain(6) -------------- elastic strains
c     stress(6) -------------- stress tensor
c     pr        -------------- Poisson's ratio
c     ym        -------------- Young's modulus
      real*8 factor,  stress(6),strain(6),ym,pr
      factor = ym/((1.+pr)*(1.-2.*pr))
      stress(1)=factor*((1.-pr)*strain(1) + pr * strain(2) + 
     $     pr * strain(3))
      stress(2)=factor*(pr*strain(1) + (1.-pr) * strain(2) + 
     $     pr * strain(3))
      stress(3)=factor*(pr*strain(1) + pr * strain(2) + 
     $     (1.-pr)*strain(3))
      stress(4)=factor*(((1.-2.*pr)/2.) * strain(4) )
      stress(5)=factor*(((1.-2.*pr)/2.) * strain(5) )
      stress(6)=factor*(((1.-2.*pr)/2.) * strain(6) )
      return
      end
c ---------------------------------------------------------VERTEX RETURN FUNCTIONS ---------------------------------------------------

      subroutine cdpm2u_checkForVertexCase(apexStress,sigV,tempkappa,
     $     rtype)
c     Subroutine that check whether the current stress state requires plasticity return to the vertex, at which
c     derivative of plastic potential and yield surface are discontinuous.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      real*8 apexStress,sigV,tempkappa,qh2,cdpm2u_qh2fun 
      integer rtype
c     sigV              -----------  volumetric stress <Input>
c     tempKappa         -----------  cummulative plastic strain  <Input>
c     apexStress        -----------  sigmaV of the yield surface for the current 
c                                    tempkappa and rho=theta=0 <Output>
c       rtype           -----------  return type of the yield surface  <Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     qh2               -----------  variables containing the results of the hardening functions 
c     cdpm2u_qh2fun     -----------  function to calculate the hardening function  given in eq. (31) of IJSS paper by P. Grassl et al.

      if ( sigV .gt. 0. ) then
         rtype = 1
         if (tempKappa .lt. 1.) then
            apexStress = 0.
         else
            qh2=cdpm2u_qh2fun(tempKappa,hp)
            apexStress=qh2*fc/m0
         end if        
      else if ( sigV .lt. 0. .and. tempKappa .lt. 1.) then
         rtype = 2
         apexStress = 0.
      else 
         rtype = 0
         apexStress=0.
      end if
      
      return 
      end
      
      
      subroutine cdpm2u_performVertexReturn(stress,apexStress,
     $     kappa,rtype,converged)
c     Subroutine that performs plasticity return close whenever the stress state needs
c     to be returned to the apex. If the stress state is not an actual vertex case
c     rtype=0 is returned.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      real*8 sigV,rho,apexStress,kappa,stress(6),theta,
     $     yieldValue,yieldValueMid,sig2,dSig,sigMid,sigAnswer,
     $     ratioPotential,kappa0,tempKappaP, cdpm2u_computeTempKappa,
     $     cdpm2u_computeRatioPotential,ratioTrial,yldTol
      integer i,j,k,rtype,maxiter,converged
c     apexStress     ------------------ variable containing the apex to return <Input>
c     kappa          ------------------ cummulative plastic strain  <Input>
c     rtype          ------------------  return type of the yield surface  <Input/Output>
c                       =0 regular return
c                       =1 return on the tensile apex of the yield surface
c                       =2 return on the compressive apex of the yield surface
c     converged      ------------------ integer denoting whether plasticity algorithm has converged
c                       =0 converged
c                       =1 not converged
c     kappa0,tempKappa---------------- variables for cummulative plastic strains
c     yieldValue,yieldValueMid ------- variables with yield values 
c     sigV           ----------------- volumetric stress 
c     rho            ----------------- deviatoric stress 
c     sig2,sigMid,sigAnswer,dSig ----- variables containing volumetric stress units
c     ratioPotential ----------------- variable containing the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     ratioTrial     ----------------- ratio of rho/sigmaV
c     yldTol         ----------------- tolerance used in bisection method solver
c     cdpm2u_computeTempKappa--------- function to calculate tempKappa according to eq.(32) of the IJSS paper by Grassl et al.
c     cdpm2u_computeRatioPotential --- function to calculate the ratio of the derivatives of plastic potential with respect to rho and sig multiplied by a parameter to convert strains in stresses
c     maxiter         ---------------- parameter denoting max number of iterations performed using the bisection method
      yldTol=gTol
      yieldValue = 0.
      yieldValueMid = 0.
      sig2 = 0.
      kappa0=kappa
      tempKappaP=kappa
      maxiter=250

      call cdpm2u_computeTrialCoordinates(stress,sigV,rho,theta)

      sig2 = apexStress
      
      tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigV,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValue,sigV, 0., 0., tempKappaP)
      
      tempKappaP =
     $     cdpm2u_computeTempKappa(kappa0, sigV, rho, sig2,ym,pr)
      
      call cdpm2u_computeYieldValue(yieldValueMid,sig2, 0.,0.,
     $     tempKappaP)
      
      if ( yieldValue * yieldValueMid .ge. 0. )  then
         converged=1
         rtype = 0  
         goto 501
      end if
      
      if ( yieldValue .lt. 0.0 ) then
         dSig = sig2 - sigV
         sigAnswer = sig2
      else 
         dSig = sigV - sig2
         sigAnswer = sig2
      end if
      
      do  j = 1, maxiter
         dSig = 0.5 * dSig
         sigMid = sigAnswer + dSig
         tempKappaP =cdpm2u_computeTempKappa(kappa0, sigV, rho, sigMid,
     $        ym,pr)
         
         call cdpm2u_computeYieldValue(yieldValueMid,sigMid, 0., 0.,
     $        tempKappaP)
         
        if ( yieldValueMid .le. 0. ) then
            sigAnswer = sigMid
         end if
         if (abs(yieldValueMid) .lt. yldTol .and. 
     $        yieldValueMid .le. 0.) then

            ratioPotential =
     $           cdpm2u_computeRatioPotential(sigAnswer, tempKappaP)
            
            ratioTrial = rho / ( sigV - sigAnswer );
            
            if ( ( ( ( ratioPotential .ge. ratioTrial ) .and. 
     $           rtype .eq. 1 ) ) .or.
     $           ( ( ratioPotential .le. ratioTrial ) .and. 
     $           rtype .eq. 2  ) ) then
               goto 500
            else    
               converged=1
               rtype = 0           
               goto 501
            endif
         endif
      enddo
 500  do k = 1, 3
         stress(k) = sigAnswer
         stress(k+3) = 0.
      enddo
      kappa=tempKappaP                          
      converged=0
 501  continue
      return
      end


      real*8 function cdpm2u_computeTempKappa(kappaInitial,sigV1,rho,
     $     sigV2,ym,pr)
c     Function to calculate the tempKappa whenever requested from the performVertexReturn function.
c     TempKappa is calculated according to eq. (32) of the IJSS paper by P. Grassl et al.

      real*8 sigV1,rho,sigV2,kappaInitial,ym,pr,
     $     equivalentDeltaPlasticStrain,kM,gM,ducMeas, 
     $     cdpm2u_computeDucMeas
c     sigV1 -------------- volumetric stress in the previous stress state <Input>
c     sigV2 -------------- volumetric stress in the current stress state <Input>
c     rho   -------------- deviatoric stress  <Input>
c     kappaInitial ------- previous kappaP (cummulative plastic strain) <Input>
c     ym    -------------- Young's modulus  <Input>
c     pr    -------------- Poisson's ratio <Input>
c     kM,gM -------------- bulk and shear moduli
c     equivalentDeltaPlasticStrain  -----  Increase of the plastic strains
c     ducMeas------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas------ function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM = ym / ( 2. * ( 1. + pr ) )

      equivalentDeltaPlasticStrain = sqrt( 1. / 9. *  (( sigV1 - sigV2 ) 
     $     /  kM )** 2.  + (rho / ( 2. * gM ))** 2. )
                                         
      ducMeas = cdpm2u_computeDucMeas(sigV2, 0., 3.141592653589793/3.)

      cdpm2u_computeTempKappa=kappaInitial+equivalentDeltaPlasticStrain/
     $     ducMeas
      return
      end

      real*8 function cdpm2u_computeRatioPotential(sig ,kappa)
c     Function to calculate the ratio of the derivatives of the plastic potential, given in eq.(22) of the IJSS paper by P. Grassl et al. with respect to the deviatoric and volumetric stress respectively.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      real*8 AGParam,BGParam,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,R,mQ,kappa,
     $     sig,rho,
     $     dgdsig,dgdrho,Al,Bl
      integer j
c     sig          ---------------- volumetric stress <Input>
c     rho          ---------------- deviatoric stress <Input>
c     kappa        ---------------- cummulative plastic strain kappaP <Input>
c     dgdsig,dgdrho --------------- derivatives of the plastic potential
c     AGParam,BGParam ------------- components of the plastic potential function
c     qh1,qh2          ------------ variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun - functions to calculate the hardening functions  given in eqs. (30) in IJSS paper by Grassl et al.
c     Al,Bl        ---------------- components of the function given in eq.(23) of the IJSS paper by P. Grassl et al.
c     R,mQ         ---------------- variables 
      rho=0.
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =    qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2.*df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(1.5) *rho/fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * (1. - qh1 ) * Bl + 6. ) +
     $     m0 * (qh1**2.) / ( sqrt(6.) * fc )
      
      cdpm2u_computeRatioPotential= dgdrho/dgdsig*3.*(1.-2.*pr ) / 
     $     ( 1. + pr )
      return 
      end

c ------------------------------------ Plasticity Algorithm functions (regular return functions) ---------------------------------
      subroutine cdpm2u_performRegularReturn(stress,kappa,converged,ym,
     $     pr,yieldTol,jacobian,normalisedResid,normOfResiduals,endflag)
c     Subroutine to perform regular plasticity return
      
      real*8 stress(6),tempKappa,ym,pr, resid(4), normalisedResid(4),PI,
     $     jacobian(4,4), inverseJac(4,4),increment(4),
     $     deltaLambda,normOfResiduals, unknowns(4),
     $     trialSig,trialRho,trialTheta,kappa,kappaP,tempKappaP,
     $     yieldTol,sig,rho,ddkappadDeltaLambdadInv(2),
     $     dgdInv(2),ddgddInv(2,2), dkappadDeltaLambda,
     $     dfdkappa, ddgdInvDkappa(2), ddkappadDeltaLambdadKappa,
     $     ddkappaddDeltaLambdadInv(2),kM,gM,dfdInv(2),sum,
     $     stressPrincipal(3),princDir(3,3),
     $     endflag
      integer i,j,iterations,totIter,converged,error
      
c     stress(6)      ------------------  effective stress components are in order xx,yy,zz,xy,yz,xz  <Input/Output>
c     kappa          ------------------  cummulative plastic strain (kappa) <Input>
c     ym             ------------------  Young's modulus  <Input>
c     pr             ------------------  Poisson's ratio <Input>
c     yieldTol       ------------------ tolerance of the N-R solver <Input>
c     converged      ------------------ integer showing whether solution has converged <Output>
c                                        = 0 solution has converged
c                                        = 1 solution has not converged
c     resid(4)       ------------------ array containing the residuals of the N-R solver
c     normalisedResid(4) -------------- array containing the normalised residuals of the N-R solver
c     PI             ------------------ parameter defined as the pi=3.14....
c     jacobian(4,4)  ------------------ matrix containing the jacobian of the problem in order to calculate the solution
c     inversJac(4,4) ------------------ matrix containing the inverse matrix of the jacobian of the problem in order to calculate the solution
c     increment(4)   ------------------ array containing the increments of the unknowns at each N-R iteration
c     deltaLambda    ------------------ plastic multiplier
c     normOfResiduals ----------------- norm the array resid(4)
c     unknowns(4)    ------------------ array with the unknowns in order sigV,rho,kappa,deltaLambda
c     trialSig       ------------------ trial volumetric stress (initial guess)
c     trialRho       ------------------ trial deviatoric stress (initial guess)
c     trialTheta     ------------------ trial Lode angle (initial guess)
c     kappaP          ------------------ cummulative plastic strain kappa (initial guess)
c     tempKappa      ------------------ temporary cummulative plastic strain kappa (each iteration)
c     sig            ------------------ temporary volumetric stress (each iteration)
c     rho            ------------------ temporary deviatoric stress (each iteration)
c     ddkappadDeltaLambdadInv(2) ------ derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric stress
c     dfdInv(2)      ------------------ derivative of the yield function with respect to the volumetric and deviatoric stress
c     dgdInv(2)      ------------------ derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     ddgddInv(2,2)  ------------------ second derivative of the plastic potential with respect to the volumetric and deviatoric stress
c     dkappadDeltaLambda -------------- derivative of kappa with respect to the plastic multiplier
c     dfdkappa       ------------------ derivative of the yield function with respect to kappa
c     ddgdInvDkappa(2) ---------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and to kappa
c     ddkappadDeltaLambdadKappa ------- derivative of kappa with respect to the plastic multiplier and kappa
c     ddkappaddDeltaLambdadInv(2) ----- derivative of kappa with respect to the plastic multiplier and deviatoric and volumetric stress
c     kM             ------------------ bulk modulus
c     gM             ------------------ shear modulus
c     sum            ------------------ variable used for summation
c     stressPrincipal(3)  ------------- array containing the principal stresses 
c     princDir(3,3)  ------------------ matrix containing eigenvectors of the effective stress tensor stored columnwise 
c     i,j,iterations ------------------ integers used as counters
c     totIter        ------------------ maximum number of iterations of the N-R algorithm
c     error          ------------------ integer indicating whether the inversion of the jacobian matrix was successful
      iterations=0
      totIter=100

      PI=3.1415926535897932384626433832795029
      kM = ym / ( 3. * ( 1. - 2. * pr ) )
      gM =  ym / ( 2. * ( 1. + pr ) )
      normOfResiduals=1.
      do i=1,4
         resid(i)=0.
         normalisedResid(i)=0.
         unknowns(i)=0.
         increment(i)=0.
      enddo
      deltaLambda=0.
      call cdpm2u_computePrincValues(stress,stressPrincipal,0,princDir)
      call cdpm2u_computeTrialCoordinates(stress,trialSig,trialRho,
     $     trialTheta)     
      kappaP=kappa
      tempKappaP=kappa
      sig=trialSig
      rho=trialRho
      unknowns(1)=trialSig
      unknowns(2)=trialRho
      unknowns(3)=tempKappaP
      unknowns(4)=0.
      
      call cdpm2u_computeYieldValue(resid(4),sig, rho,trialTheta,
     $     tempKappaP)
      normOfResiduals=1.
      do while (normOfResiduals .gt. yieldTol)
c      write(*,*) 'normOfResiduals = ',normOfResiduals
         iterations=iterations+1
         if (iterations .eq. totIter) then
            converged=1
            goto 600
         end if 
         normalisedResid(1)=resid(1)/kM
         normalisedResid(2)=resid(2)/2./gM
         normalisedResid(3)=resid(3)
         normalisedResid(4)=resid(4)
         normOfResiduals=sqrt(normalisedResid(1)**2.+
     $        normalisedResid(2)**2.+normalisedResid(3)**2. +
     $        normalisedResid(4)**2.)
         
         if (isnan(normOfResiduals)) then
            converged=1
            goto 600
         end if

         if (normOfResiduals .gt. yieldTol) then
c     ----------------- compute jacobian ---------------------------------
           call cdpm2u_computedfdInv(dfdInv,sig,rho,trialTheta,
     $           tempKappaP) 
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgddInv(ddgddInv,sig,rho,trialTheta,
     $          tempKappaP)
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP,endflag)
           call cdpm2u_computedfdKappa(dfdkappa,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $          trialTheta,tempKappaP,endflag)
           call cdpm2u_computeddKappadDeltaLambdadKappa(
     $          ddkappadDeltaLambdadKappa,sig,rho,tempKappaP,trialTheta,
     $          endflag)

           call cdpm2u_computeddKappadDeltaLambdadInv(
     $          ddKappaddDeltaLambdadInv,sig,rho,tempKappaP,trialTheta,
     $          endflag)

           jacobian(1,1) = 1. + kM * deltaLambda *   ddgddInv(1, 1)
           jacobian(1, 2) = kM * deltaLambda * ddgddInv(1, 2)
           jacobian(1, 3) = kM * deltaLambda * ddgdInvdKappa(1)
           jacobian(1, 4) = kM * dgdInv(1)
           
           jacobian(2, 1) = 2. *gM *deltaLambda *ddgddInv(2, 1)
           jacobian(2, 2) = 1. + 2. *gM *deltaLambda *  ddgddInv(2, 2)
           jacobian(2, 3) = 2. *gM *deltaLambda * ddgdInvdKappa(2)
           jacobian(2, 4) = 2. *gM *dgdInv(2)
           
           jacobian(3, 1) = deltaLambda * ddKappaddDeltaLambdadInv(1)
           jacobian(3, 2) = deltaLambda * ddKappaddDeltaLambdadInv(2)
           jacobian(3, 3) = deltaLambda * ddkappadDeltaLambdadKappa - 1.
           jacobian(3, 4) = dkappadDeltaLambda
           
           jacobian(4, 1) = dfdInv(1)
           jacobian(4, 2) = dfdInv(2)
           jacobian(4, 3) = dfdKappa
           jacobian(4, 4) = 0.

           call cdpm2u_computeInverseJac(inverseJac,jacobian,error)
           if (error.eq.-1) then
              converged=1
              goto 600
           end if
           do i=1,4
              sum = 0.
              do j = 1,4
                 sum =sum+ inverseJac(i, j) * resid(j)
              enddo
              increment(i) =0.-sum
              unknowns(i)=unknowns(i)+increment(i)
           enddo
           if (unknowns(4) .le. 0.) then
              unknowns(4)=0.
           end if
           if (unknowns(2) .le. 0.) then
              unknowns(2)=0.
           end if
           if (unknowns(3)-kappaP .le. 0.) then
              unknowns(3)=kappaP
           end if
           sig = unknowns(1)
           rho = unknowns(2)
           tempKappaP=unknowns(3)
           deltaLambda=unknowns(4)
           call cdpm2u_computedgdInv(dgdInv,sig,rho,trialTheta,
     $          tempKappaP) 
           call cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,sig,
     $          rho,trialTheta,tempKappaP,endflag)
           resid(1) = sig - trialSig + kM *deltaLambda * dgdInv(1)
           resid(2) = rho - trialRho +  2.* gM *deltaLambda * dgdInv(2)
           resid(3) = -tempKappaP +kappaP+deltaLambda*dkappadDeltaLambda
           call cdpm2u_computeYieldValue(resid(4),sig,rho,trialTheta,
     $          tempKappaP)
        end if
      end do
      converged=0
      
      
      stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(trialTheta)
      stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta - 2. * PI/ 3.)
      stressPrincipal(3) = sig + sqrt(2. / 3.) * rho * 
     $     cos(trialTheta + 2. * PI / 3.)

      call cdpm2u_transformStressVectorTo(stress,princDir,
     $     stressPrincipal)
      
      kappa=tempKappaP
 600  continue
 
c      write(*,*) 'resid = ',resid
c      write(*,*) 'normOfResiduals = ',normOfResiduals
c      write(*,*) 'jacobian(:,1)', jacobian(:,1)
c      write(*,*) 'jacobian(:,2)', jacobian(:,2)
c      write(*,*) 'jacobian(:,3)', jacobian(:,3)
c      write(*,*) 'jacobian(:,4)', jacobian(:,4)
c      write(*,*) '----------------------------------------'
 
      return
      end
      
c ------------------------------------ Plasticity Algorithm functions (general functions) ----------------------------------------
      subroutine cdpm2u_computedfdInv(dfdInv,sig,rho,theta,kappa) 
c     Subroutine to calculate the derivative of the yield function with respect to the volumetric and deviatoric stresses respectively 
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 rFunction,qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     dfdsig,dfdrho,sig,dfdInv(2),Bl
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdInv(2)        -------------  derivatives of yield function with respect to volumeric and deviatoric stress <Output>
c     dfdsig           -------------  derivative of the yield function with respect to volumetric stress
c     dfdrho           -------------  derivative of the yield function with respect to deviatoric stress
c     Al,Bl            -------------  variables corresponding to components of the yield function
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.       
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      
      
      Al =( 1. - qh1 ) * ( sig / fc + rho / ( sqrt(6.) * fc ) )** 2. +
     $     sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      
      dfdsig= 4. * ( 1. - qh1 ) / fc * Al * Bl + qh2* qh1**2. * m0 / fc
      dfdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6. )+ 
     $     rFunction * m0 * qh2 * qh1**2./ ( sqrt(6.) * fc )
      
      dfdInv(1) = dfdsig
      dfdInv(2) = dfdrho
      return
      end
      
      subroutine cdpm2u_computedgdInv(dgdInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the plastic potential function with respect to volumetric and deviatoric stress
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,cdpm2u_qh2fun,kappa,
     $     Bl,AGParam,BGParam,R,mQ,dgdsig,dgdrho,sig,dgdInv(2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2)        -------------  derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     dgdsig           -------------  derivative of the plastic potential with respect to volumetric stress
c     dgdrho           -------------  derivative of the plastic potential with respect to deviatoric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.      
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) /  ( log(AGParam)+
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      
      dgdsig = 4. * ( 1. - qh1 ) / fc * Al * Bl + qh1**2. * mQ / fc
      dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 ) * Bl + 6.) + 
     $     m0 * qh1**2. / ( sqrt(6.) * fc )
      
      dgdInv(1) = dgdsig
      dgdInv(2) = dgdrho
      return 
      end 
      
      subroutine cdpm2u_computeddgddInv(ddgddInv,sig,rho,theta,kappa)
c     Subroutine to calculate the derivatives of the derivatives of the plastic potential with respect to volumetric and deviatoric stress
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 qh1,qh2, Al,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa, Bl,AGParam,
     $     BGParam,R,mQ,sig, dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho,
     $     ddgddSig,ddgddRho,ddgdSigdRho,ddgdRhodSig,ddgddInv(2,2)
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dgdInv(2,2)      -------------  derivatives of derivatives of plastic potential with respect to volumeric and deviatoric stress <Output>
c     ddgddsig         -------------  second derivative of the plastic potential with respect to volumetric stress
c     ddgddrho         -------------  second derivative of the plastic potential with respect to deviatoric stress
c     ddgdsigdrho      -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     ddgdrhodsig      -------------  derivative of the plastic potential with respect to deviatoric and volumetric stress
c     Al,Bl,R,mQ,AGParam,BGParam ---  variables corresponding to components of the plastic potential
c     dMQDSig,dAlDSig,dBlDSig,dAlDRho, dBlDRho, ---  variables corresponding to derivatives of components of the plastic potential
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun --  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.            
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      
      AGParam = ft * qh2 * 3. / fc + m0 / 2.
      BGParam =qh2 / 3. * ( 1. + ft / fc ) / ( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2 + m0 / 2.) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * Bl**2. + sqrt(3. / 2.) * rho / fc
      dMQDSig = AGParam / ( BGParam * fc ) * exp(R)
      dAlDSig = 2. * ( 1. - qh1 ) * Bl / fc
      dBlDSig = 1. / fc
      dAlDRho = 2. * ( 1. - qh1 ) * Bl / ( fc * sqrt(6.) ) + 
     $     sqrt(3. / 2.) / fc;
      dBlDRho = 1. / ( fc * sqrt(6.) )
      
      ddgddSig = 4. * ( 1. - qh1 ) / fc * ( dAlDSig * Bl + Al * 
     $     dBlDSig ) +qh1**2. * dMQDSig / fc
      ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * Bl + 6. ) +Al * dBlDRho * 4. *
     $     ( 1. - qh1 ) / ( sqrt(6.) * fc )
      ddgdSigdRho = 4. * (1. - qh1 )/fc *( dAlDRho * Bl + Al * dBlDRho )
      ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - 
     $     qh1 ) * Bl + 6. ) + Al / ( sqrt(6.) * fc ) * ( 4. * 
     $     ( 1. - qh1 ) * dBlDSig )
      
      ddgddInv(1, 1) = ddgddSig
      ddgddInv(1, 2) = ddgdSigdRho
      ddgddInv(2, 1) = ddgdRhodSig
      ddgddInv(2, 2) = ddgddRho
      
      return
      end
      
      subroutine cdpm2u_computedkappadDeltaLambda(dkappadDeltaLambda,
     $     sig,rho,theta,kappa,endflag)
c     Subroutine to calculate the derivative of the hardening variable with respect to plastic multiplier
      real*8 rho,sig,theta,kappa,dkappadDeltaLambda,equivalentDGDStress,
     $     ductilityMeasure,dgdInv(2),cdpm2u_computeDucMeas,endflag
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     dkappadDeltaLambda -----------  derivative of the hardening variable with respect to plastic multiplier <Output>
c     dgdInv(2)        -------------  derivative of plastic potential with respect to volumeric and deviatoric stress
c     ductilityMeasure -------------  ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- norm of the derivative of the plastic potential with respect to volumetric and deviatoric stress
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      equivalentDGDStress = sqrt( 1. / 3.*dgDInv(1)** 2.+dgdInv(2)** 2.)
      ductilityMeasure = cdpm2u_computeDucMeas(sig, rho,theta)
	  dkappadDeltaLambda = equivalentDGDStress / ductilityMeasure
	  
      if (endflag .eq. 1) then
         write(*,*) '--------------------------------------------------------'
         write(*,*) 'dkappadDeltaLambda', dkappadDeltaLambda
         write(*,*) 'sqrt( 1. / 3.*dgDInv(1)** 2.+dgdInv(2)** 2.) / ductilityMeasure', sqrt( 1. / 3.*dgDInv(1)** 2.+dgdInv(2)** 2.) / ductilityMeasure
         write(*,*) 'equivalentDGDStress', equivalentDGDStress
         write(*,*) 'dgDInv', dgDInv
         write(*,*) 'ductilityMeasure', ductilityMeasure
         write(*,*) '--------------------------------------------------------'
      end if
	  
      return
      end
      
      subroutine cdpm2u_computedfdKappa(dfdkappa,sig,rho,theta,kappa)
c     Subroutine to calculate the derivative of the yield function with respect to the cummulative plastic strain(kappa)
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      real*8 dfdkappa,sig,rho,theta,kappa,qh1,qh2,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,dfdqh1,
     $     dfdqh2,dqh1dkappa,dqh2dkappa,cdpm2u_dqh1dkappaFun,
     $     cdpm2u_dqh2dkappaFun,Al,
     $     Bl,rFunction
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     dfdkappa         -------------  derivative of the yield function with respect to cummulative plastic strain (kappa) <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     Al,Bl            -------------  variables corresponding to components of the yield function
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *cos(theta)**2.
     $     + 5. * ecc**2. - 4. * ecc) )
      Al = ( 1. - qh1 ) * ( ( sig / fc + rho / ( sqrt(6.) * 
     $     fc ) )) **2.  + sqrt(3. / 2.) * rho / fc
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      dfdqh1 = -2. *Al *(Bl** 2.) + 2. * qh1 * qh2 *   m0 * ( sig / fc +
     $     rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *qh1*(qh2**2.)

      dfdqh2 = (qh1**2.) * m0 * ( sig / fc + rho * rFunction /
     $     (sqrt(6.) * fc)) -  2. *qh2 *(qh1** 2.)
      dfdkappa =  dqh1dkappa * dfdqh1 + dqh2dkappa * dfdqh2
      
      if ( dfdkappa .gt. 0. ) then
         dfdkappa = 0.
      end if
      
      return
      end      
      
      subroutine cdpm2u_computeddgdInvdKappa(ddgdInvdKappa,sig,rho,
     $     theta,kappa,endflag)
c     Subroutine to calculate the derivative of the plastic potential function with respect to the volumetric and deviatoric stresses and the cummulative plastic strain (kappa)      
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag,
     2     endflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 qh1,cdpm2u_qh1fun,qh2,cdpm2u_qh2fun,dqh1dkappa,
     $     cdpm2u_dqh1dkappaFun,dqh2dkappa,
     $     cdpm2u_dqh2dkappaFun, AGParam,BGParam,R,mQ,dAGParamdKappa,
     $     BGParamTop,BGParamBottom,dBGParamTopDKappa,
     $     dBGParamBottomDKappa,dBGParamDKappa,RTop,RBottom,dRTopDKappa,
     $     dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard,
     $     dDGDSigDKappa,ddgdInvdKappa(2),kappa,sig,rho,theta,
     $     dDGDRhoDKappa
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain        <Input>
c     ddgdInvdKappa(2) -------------  derivative of the plastic potential with respect to deviatoric and volumetric strains and the cummulative plastic strain (kappa) <Output>
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.     
c     dqh1dkappa,dqh2dkappa --------  variables containing the results of the derivatives of the hardening functions with respect to cummulative plastic strain (kappa) 
c     cdpm2u_dqh1dkappaFun,cdpm2u_dqh2dkappaFun --  functions to calculate the derivatives of the hardening functions with respect to the  cummulative plastic strain (kappa) 
c     dDGDSigDKappa    -------------  derivative of the plastic potential with respect to volumetric stress and cummulative plastic strain (kappa)
c     dDGDRhoDKappa    -------------  derivative of the plastic potential with respect to deviatoric stress and cummulative plastic strain (kappa)
c     Al,Bl,AGParam,BGParam,R,mQ,dAGParamdKappa,BGParamTop,BGParamBottom,dBGParamTopDKappa,  dRBottomDKappa,dRDKappa,dMQDKappa,Al,Bl,dAlDYieldHard           -------------  variables corresponding to components and their derivatives of the plastic potential
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)
      dqh1dkappa=cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
      dqh2dkappa=cdpm2u_dqh2dkappaFun(kappa,hp)
      
      AGParam = ft * qh2 * 3 / fc + m0 / 2
      BGParam =  qh2 / 3. * ( 1. + ft / fc ) /( log(AGParam) + 
     $     log(df + 1.) - log(2 * df - 1.) - log(3. * qh2+ m0 / 2) )
      R = ( sig - ft / 3. * qh2 ) / fc / BGParam
      mQ = AGParam * exp(R)
      dAGParamDKappa = dqh2dkappa * 3. * ft / fc
      BGParamTop = qh2 / 3. * ( 1. + ft / fc );
      BGParamBottom = ( log(AGParam) + log(df + 1.) - 
     $     log(2 * df - 1.) - log(3. * qh2 + m0 / 2) )
      dBGParamTopDKappa = dqh2dkappa / 3.
      dBGParamBottomDKappa = -3. * dqh2dkappa / ( 3 * qh2 + m0 / 2. )
      dBGParamDKappa =( dBGParamTopDKappa * BGParamBottom - BGParamTop * 
     $     dBGParamBottomDKappa ) / (BGParamBottom**2.)
      RTop = ( sig - ft / 3. * qh2 )
      RBottom = fc * BGParam
      dRTopDKappa = -ft / 3. * dqh2dkappa
      dRBottomDKappa = fc * dBGParamDKappa
      dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / 
     $     (RBottom** 2.)
      dMQDKappa = dAGParamDKappa * exp(R) + AGParam *dRDKappa *exp(R)
      Bl = sig / fc + rho / ( fc * sqrt(6.) )
      Al = ( 1. - qh1 ) * (Bl** 2.) + sqrt(3. / 2.) * rho / fc
      dAlDYieldHard = -(Bl** 2.)
      
      dDGDSigDKappa =  ( -4. * Al * Bl / fc + 4. * ( 1 - qh1 ) / fc * 
     $     dAlDYieldHard * Bl ) * dqh1dkappa +
     $     dqh1dkappa * 2 * qh1 * mQ / fc + qh1 * dMQDKappa / fc
      dDGDRhoDKappa =
     $     ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - qh1 )* 
     $     Bl + 6. ) - 4. * Al / ( sqrt(6.) * fc ) * Bl + m0/( sqrt(6.)* 
     $     fc ) ) * 2 * qh1 * dqh1dkappa
      
      ddgdInvdKappa(1) = dDGDSigDKappa
      ddgdInvdKappa(2) = dDGDRhoDKappa
	  
      if (endflag .eq. 1) then
         write(*,*) '--------------------------------------------------------'
         write(*,*) 'ddgdInvdKappa', ddgdInvdKappa
         write(*,*) 'dAlDYieldHard', dAlDYieldHard
         write(*,*) 'Al', Al
         write(*,*) 'Bl', Bl
         write(*,*) 'qh1', qh1
         write(*,*) 'dqh1dkappa', dqh1dkappa
         write(*,*) 'kappa', kappa
         write(*,*) '--------------------------------------------------------'
      end if
      
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadKappa(
     $     ddkappadDeltaLambdadKappa,sig,rho,kappa,theta,endflag)
c     Subroutine to compute the derivative of the cummulative plastic strain(kappa) with respect to the plastic multiplier and the hardening variable kappa
      real*8 equivalentDGDStress,dEquivalentDGDStressDKappa,ducMeas,
     $     ddkappadDeltaLambdadKappa,dgdInv(2),ddgdInvdKappa(2),
     $     sig,rho,kappa,theta,cdpm2u_computeDucMeas,endflag,
     $     ddkappadDeltaLambdadKappa2
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddkappadDeltaLambdadKappa ----  derivative of the kappa with respect to plastic multiplier and kappa <Output> 
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     ddgdInvdKappa(2) ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress and cummulative plastic strain (kappa)
c     dgdInv(2)        ------------- derivative of the plastic potential with respect to volumetric and deviatoric stress
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDKappa ---  scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress and kappa which is the norm of the derivative of the increment of the plastic strains with respect to kappa 
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgdInvdKappa(ddgdInvdKappa, sig, rho,theta,
     $     kappa,endflag)      
      equivalentDGDStress =sqrt( 1./ 3.*( dgdInv(1)**2.)+ dGDInv(2)**2.)
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivalentDGDStressDKappa = (2./3.*dgdInv(1)/equivalentDGDStress*
     $     ddgdInvdKappa(1) + 2.*dgdInv(2)/equivalentDGDStress*
     $     ddgdInvdKappa(2)) / 2.
      ddkappadDeltaLambdadKappa= dEquivalentDGDStressDKappa/ducMeas
	  
      if (endflag .eq. 1) then
         write(*,*) '--------------------------------------------------------'
         write(*,*) 'ddkappadDeltaLambdadKappa', ddkappadDeltaLambdadKappa
         write(*,*) 'ddkappadDeltaLambdadKappa2', ddkappadDeltaLambdadKappa2
         write(*,*) 'dEquivalentDGDStressDKappa/ducMeas', dEquivalentDGDStressDKappa/ducMeas
         write(*,*) 'equivalentDGDStress', equivalentDGDStress
         write(*,*) 'dgdInv', dgdInv
         write(*,*) 'dEquivalentDGDStressDKappa', dEquivalentDGDStressDKappa
         write(*,*) 'ddgdInvdKappa', ddgdInvdKappa
         write(*,*) 'ducMeas', ducMeas
         write(*,*) '--------------------------------------------------------'
      end if
      return
      end
      
      subroutine cdpm2u_computeddKappadDeltaLambdadInv(
     $     ddKappaddDeltaLambdadInv, sig,rho,kappa,theta,endflag)
c     Subroutine to compute the derivative of the cummulative plastic strain (kappa) with respect to the plastic multiplier and the volumetric and deviatoric stress
      real*8 equivDGDStress, dgdInv(2),ddgddInv(2, 2),
     $     dEquivDGDStressDInv(2),ducMeas,
     $     ddKappaddDeltaLambdadInv(2),sig,rho,kappa,theta,
     $     dDuctilityMeasureDInv(2),cdpm2u_computeDucMeas,endflag
c     sig              -------------  volumetric stress  <Input>
c     rho              -------------  deviatoric stress  <Input>
c     theta            -------------  Lode  angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     ddKappaddDeltaLambdadInv -----  derivative of the kappa with respect to plastic multiplier and volumetric and deviatoric strain <Output> 
c     ddgdd(2,2)       -------------  second derivative of the plastic potential with respect to volumetric and deviatoric stress 
c     dgdInv(2)        -------------  derivative of the plastic potential with respect to volumetric and deviatoric stress
c     dDuctilityMeasureDInv(2) -----  derivative of the ductility measure with respect to volumetric and deviatoric stress
c     ducMeas          ------------- ductility measure in plasticity
c     cdpm2u_computeDucMeas   ------------- function to calculate the ductility measure according to eq.(33) of IJSS paper by P. Grassl et al.
c     equivalentDGDStress ---------- scalar function of the derivative of the plastic potential with respect to volumetric and deviatoric stress which is the norm of the increment of the plastic strains 
c     dEquivalentDGDStressDInv(2) -- scalar function of the second derivative of the plastic potential with respect to volumetric and deviatoric stress which is the derivative of the norm of the increment of the plastic strains with respect to deviatoric and volumetric plastic strains 
      
      call cdpm2u_computedgdInv(dgdInv, sig, rho,theta, kappa)
      call cdpm2u_computeddgddInv(ddgddInv, sig, rho,theta, kappa)
      
      equivDGDStress = sqrt( 1. / 3. * (dGDInv(1)** 2.) +
     $     dGDInv(2)**2. )
      ducMeas = cdpm2u_computeDucMeas(sig, rho, theta)
      dEquivDGDStressDInv(1) =  ( 2. / 3.*dgdInv(1)*ddgddInv(1,1) + 
     $     2. * dgdInv(2)*ddgddInv(2, 1) ) / ( 2. * equivDGDStress)
      dEquivDGDStressDInv(2) = ( 2. / 3.*dgdInv(1)*ddgddInv(1,2) +
     $     2. * dgdInv(2)*ddgddInv(2, 2) ) /( 2. * equivDGDStress )
      call cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, rho,
     $     theta, kappa)
      
      ddKappaddDeltaLambdadInv(1) = ( dEquivDGDStressDInv(1) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(1) ) / (ducMeas** 2.)
      ddKappaddDeltaLambdadInv(2) = ( dEquivDGDStressDInv(2) * ducMeas - 
     $     equivDGDStress * dDuctilityMeasureDInv(2) ) / (ducMeas** 2.)	 
      
      if (endflag .eq. 1) then
         write(*,*) '--------------------------------------------------------'
         write(*,*) 'ddKappaddDeltaLambdadInv', ddKappaddDeltaLambdadInv
         write(*,*) 'dEquivDGDStressDInv(1) * ducMeas', dEquivDGDStressDInv(1) * ducMeas
         write(*,*) 'equivDGDStress * dDuctilityMeasureDInv(1)', equivDGDStress * dDuctilityMeasureDInv(1)
         write(*,*) 'dEquivDGDStressDInv(2) * ducMeas', dEquivDGDStressDInv(2) * ducMeas
         write(*,*) 'equivDGDStress * dDuctilityMeasureDInv(2)', equivDGDStress * dDuctilityMeasureDInv(2)
         write(*,*) 'dEquivDGDStressDInv', dEquivDGDStressDInv
         write(*,*) 'ducMeas', ducMeas
         write(*,*) 'equivDGDStress', equivDGDStress
         write(*,*) 'dDuctilityMeasureDInv', dDuctilityMeasureDInv
         write(*,*) 'dgdInv', dgdInv
         write(*,*) 'ddgddInv', ddgddInv
         write(*,*) '--------------------------------------------------------'
      end if	  
      return
      end
      
      
      subroutine cdpm2u_computeTrialCoordinates(stress, sigV,rho, theta)
c     Subroutine which returns volumetric and deviatoric stress and the Lode angle
c     based on the given stress tensor
      real*8 stress(6),rho,sigV,theta,tempdevSig(6),j2,j3
      integer i,j
      
c     stress(6)     -------------  Stress array. stress={sigXX,sigYY,sigZZ,sigXY,sigYZ,sigXZ} <Input>
c     sigV          -------------  Volumetric stress. <Output>
c     rho           -------------  Norm of deviatoric stress. <Output>
c     theta         -------------  Lode angle.  <Output>
c     tempdevSig(6) -------------  Array containing  deviatoric stress tensor
c     j2            -------------  Second invariant of the deviatoric stress tensor
c     j3            -------------  Third invariant of the deviatoric stress tensor
c     i,j           -------------  counters used in iterations
      
      sigV= (stress(1)+stress(2)+stress(3))/3.
      do j=1,3
         tempdevSig(j)=stress(j)-sigV
         tempdevSig(j+3)=stress(j+3)
      enddo
      j2=0.5*(tempdevSig(1)*tempdevSig(1)+tempdevSig(2)*tempdevSig(2)+ 
     $     tempdevSig(3)*tempdevSig(3))+ tempdevSig(4)*tempdevSig(4)+ 
     $     tempdevSig(5)*tempdevSig(5)+tempdevSig(6)*tempdevSig(6)
      
      if(j2 .eq. 0.) then
         theta=0.
         rho=0.
      else   
         rho=sqrt(2.*j2)
         j3= (1./3.) * ( tempdevSig(1)*tempdevSig(1)*tempdevSig(1) + 
     $        3.*tempdevSig(1) *tempdevSig(4)*tempdevSig(4)+
     $        3.*tempdevSig(1) * tempdevSig(6) * tempdevSig(6) + 
     $        6. * tempdevSig(5) * tempdevSig(4)  *  tempdevSig(6) +
     $        3. * tempdevSig(2) * tempdevSig(4)**2.+
     $        3 * tempdevSig(3) *tempdevSig(6)*tempdevSig(6)+
     $        tempdevSig(2)*tempdevSig(2)*tempdevSig(2)+ 
     $        3. * tempdevSig(2) * tempdevSig(5)*tempdevSig(5)+
     $        3. * tempdevSig(3)*tempdevSig(5)*tempdevSig(5)+ 
     $        tempdevSig(3)*tempdevSig(3)*tempdevSig(3))
         theta=(3.*sqrt(3.)/2.)*j3/(j2**(3./2.))
      end if
      if (theta .gt. 1.) then
         theta=1.
      else  if (theta .lt. -1.) then
         theta=-1.
      end if
      theta=1./3.*acos(theta)
      return
      end
      
      subroutine cdpm2u_computeYieldValue(answer,sigV, rho,theta,kappa)
c     Function to evaluate the yield function f based on the given stress High -Westergaard coordinates and kappa.
c     The equation is given in eq. (18) of IJSS paper by P. Grassl et al. Returns a real*8 number 
      
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      
      real*8 rFunction,qh1,qh2, Al,sigV,rho,theta,cdpm2u_qh1fun,
     $     cdpm2u_qh2fun,kappa,
     $     answer

c     sigV             -------------  volumetric stress <Input>
c     rho              -------------  deviatoric stress <Input>
c     theta            -------------  Lode angle        <Input>
c     kappa            -------------  cummulative plastic strain  <Input>
c     answer           -------------  variable containing the result of the yield function <Output>
c     rFunction        -------------  function to control shape of the yield surface given in eq. (19) of IJSS paper by P. Grassl et al.
c     qh1,qh2          -------------  variables containing the results of the hardening functions 
c     cdpm2u_qh1fun,cdpm2u_qh2fun    -------------  functions to calculate the hardening functions  given in eqs. (30), (31) of IJSS paper by P. Grassl et al.
c     Al               -------------  variable used to siplify and facilitate the calculation of the yield surface
      qh1=cdpm2u_qh1fun(kappa,qh0,hp)
      qh2=cdpm2u_qh2fun(kappa,hp)

      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(theta)**2. +
     1     ( 2. * ecc - 1. )**2.  ) /
     2     ( 2. * ( 1. - ecc**2. ) * cos(theta) +
     3     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     4     cos(theta)**2.+ 5. * ecc**2. - 4. * ecc) )                                                              

      Al = ( 1. - qh1 ) * (sigV / fc + rho / (sqrt(6.) * fc) ) ** 2. +
     &     sqrt(1.5) * rho / fc
      
      answer=Al**2. + qh1**2. * qh2* m0 * ( sigV / fc + rho *
     &     rFunction / ( sqrt(6.) * fc ) ) - qh1**2. * qh2**2.      
      return
      end
      
      
      real*8 function cdpm2u_qh1fun(kappa,qh0,hp)
c     Function to calculate the first hardening function given in eq. (30) of the IJSS paper 
c     by Grassl et al.

      real*8 hp,qh0, kappa,answer

c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     qh0    ---------- Initial hardening parameter <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer=qh0
      else if ((kappa .gt. 0.0) .and. (kappa .lt. 1.0)) then
         answer=
     $        (1.0 - qh0 - hp) *(kappa**3.0)-
     $        ( 3.0 * ( 1.0 - qh0 ) - 3.0 * hp ) * (kappa**2.0) +
     $        ( 3.0 * ( 1.0 - qh0 ) - 2.0 * hp ) * kappa + qh0
      else 
         answer=1.
      endif      
      cdpm2u_qh1fun=answer      
      return
      end
      
      real*8 function cdpm2u_qh2fun(kappa,hp)
c     Function to calculate the first hardening function given in eq. (31) of the IJSS paper 
c     by Grassl et al.
      
      real*8 kappa,answer,hp
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      
      if ( kappa .le. 0. ) then
         answer=1.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=1.
      else 
         answer= 1.+(kappa-1.)*hp
      endif      
      cdpm2u_qh2fun=answer     
      return 
      end
      
      real*8 function cdpm2u_dqh1dkappaFun(kappa,qh0,hp)
c     Function to calculate the derivative of the first hardening function, given in eq. (30) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real*8 kappa,answer,hp,qh0
c     kappa  ---------- Cummulative plastic strain <Input> 
c     qh0    ---------- Initial hardening parameter <Input>
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function

      if ( kappa .le. 0. ) then
         answer= 3. * ( 1 - qh0 ) - 2. * hp
      else if ( kappa .ge. 0. .and. kappa .lt. 1. ) then
         answer=   3. * ( 1. - qh0 - hp ) * (kappa**2.)
     $        - 2. * ( 3. * ( 1. - qh0 ) - 3. * hp ) * kappa
     $        + ( 3. * ( 1. - qh0 ) - 2. * hp )
      else 
         answer=  0.
      endif
      cdpm2u_dqh1dkappaFun=answer
      return
      end
      
      real*8 function cdpm2u_dqh2dkappaFun(kappa,hp)
c     Function to calculate the derivative of the second hardening function, given in eq. (31) of the IJSS paper 
c     by Grassl et al., with respect to the cummulative plastic strain (kappaP)
      real*8 kappa,hp,answer
c     kappa  ---------- Cummulative plastic strain <Input> 
c     hp     ---------- Hardening modulus <Input>
c     answer ---------- Variable containing the answer of the function
      if ( kappa .le. 0. ) then
         answer=0.
      else if ( kappa .gt. 0. .and. kappa .lt. 1. ) then
         answer=0.
      else 
         answer=hp
      endif
      cdpm2u_dqh2dkappaFun=answer
      return 
      end  


      real*8 function cdpm2u_computeDucMeas( sigV,rho,theta)
c     Function to calculate ductility measure according to eq. (33) of IJSS paper by P. Grassl et al.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      real*8 thetaConst,x,  sigV,rho,theta ,eh,fh,answer
c     sigV ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta---------------- Lode angle <Input>
c     eh,fh---------------- hardening parameters
c     thetaConst,x -------- variables
c     answer--------------- variable containing the answer of the function

      thetaConst = (2. * cos(theta))**2.
      x = -( sigV + fc / 3 ) / fc
      if ( x .lt. 0. ) then
         eh = bh - dh
         fh = ( bh - dh ) * ch / ( ah - bh )
         answer = ( eh * exp(x / fh) + dh ) / thetaConst
      else 
         answer = ( ah + ( bh - ah ) * exp( -x / ( ch ) ) ) / thetaConst
      endif
      cdpm2u_computeDucMeas=answer
      return
      end
      
      
      
      subroutine cdpm2u_computedDucMeasdInv(dDuctilityMeasureDInv, sig, 
     $     rho,theta, kappa)
c     Subroutine to compute the derivative of the ductility measure with respect to volumetric and deviatoric stress
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 dDuctilityMeasureDInv(2), sig,rho, kappa,theta,x,dXDSig,
     $     EHard,FHard,dDuctilityMeasureDX,theta1
c     sig  ---------------- volumetric stress <Input>
c     rho  ---------------- deviatoric stress <Input>
c     theta --------------- Lode angle <Input>
c     dDuctilityMeasureDInv(2) --  derivative of the ductility measure with respect to volumetric and deviatoric stress <Output>
c     eh,fh --------------- hardening parameters
c     theta1,x ------------ variables
c     EHard,FHard --------- hardening parameters
c     dXDSig -------------- derivative of the function R, given in eq. (34) in IJSS paper by P.Grassl et al., with respect to volumetric stress
c     dDuctilityMeasureDX - derivative of the ductility measure with respect to the function R

      theta1 = (2. * cos(theta))** 2.
      x = ( -( sig + fc / 3. ) ) / fc
      
      if ( x .lt. 0. ) then
         dXDSig = -1. / fc
         EHard = bh - dh
         FHard = ( bh - dh ) * ch / ( ah - bh )
         
         dDuctilityMeasureDX = EHard  / FHard *exp(x / FHard) / theta1
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      else 
         dXDSig = -1. / fc
         dDuctilityMeasureDX = -( bh - ah ) / ( ch ) / theta1 *
     $        exp( -x / ( ch ) )
         dDuctilityMeasureDInv(1) = dDuctilityMeasureDX * dXDSig
         dDuctilityMeasureDInv(2) = 0.
      endif
      return
      end
c     --------------------------------------------- Functions used in the damage algorithm ---------
      subroutine cdpm2u_computeAlpha(stressT,stressC,stress,alpha)
      
c     Subroutine calculating the alpha according to eq.46 of the IJSS paper by Grassl et al. and splitting stress tensor in tensile and compressive part in the principal effective stress coordinate system. Finally principal effective tensile and compressive stress tensors are rotated back to the original system

      real*8 stress(6),stressT(6),stressC(6)
      real*8 princDir(3, 3),princStress(3),princStressT(3),alpha,
     $     princStressC(3), squareNormOfPrincipalStress,alphaTension
      integer i
c     stress(6)       -------------  effective stress tensor (no damage included) <Input>
c     stressT(6)      -------------  effective tensile stress tensor (no damage included) <Output>
c     alpha           ------------- alphaC used in the material model <Output>
c     stressC(6)      -------------  effective compressive stress tensor (no damage included) 
c     alphaTension    ------------- 1-alphaC <Output>
c     princDir(3,3)   -------------  matrix containing the principal directions of original stress tensor stored columnwise
c     princStress(3)  ------------- array containing principal stresses {sigma1,sigma2,sigma3} 
c     princStressT(3) ------------- array containing principal tensile stresses {sigmaT1,sigmaT2,sigmaT3} 
c     princStressC(3) ------------- array containing principal compressive stresses {sigmaC1,sigmaC2,sigmaC3} 
c     squareNormOfPrincipalStress - norm of princStress 

      call cdpm2u_computePrincValues(stress,princStress,0,princDir)
      
c     Split the principal values in a tension and a compression part
      do i = 1,3
         if ( princStress(i) .ge. 0. ) then
            princStressT(i) = princStress(i)
            princStressC(i) = 0.
         else 
            princStressC(i) = princStress(i)
            princStressT(i) = 0.
         endif
      enddo
      
c     Transform the tension and compression principal stresses back to the original coordinate system

      call cdpm2u_transformStressVectorTo(stressT, princDir,
     $     princStressT)
      call cdpm2u_transformStressVectorTo(stressC, princDir,
     $     princStressC)
      
c     Determine the two factors from the stress
      squareNormOfPrincipalStress = 0.

         squareNormOfPrincipalStress = princStress(1)**2.+
     $     princStress(2)**2.+princStress(3)**2.

      
      alphaTension = 0.
      
      if ( squareNormOfPrincipalStress .gt. 0. ) then
         do i = 1,3
            alphaTension = alphaTension+ princStressT(i) *
     $           ( princStressT(i) + princStressC(i) ) /
     $           squareNormOfPrincipalStress
         enddo
      endif

      alpha= 1. - alphaTension

      return
      end

      real*8 function cdpm2u_computeRateFactor(alpha,strainrate,fc,
     $     fc0,tol)
c     Function to incorporate impact effects in the constitutive law and calculate rateFactor. 
c     All functions used are based on the equations of the chapter 2.1.5 of the Model Code 1990. 
      real*8 strainrate(6),princDir(3,3),princStrainRate(3),max,min,tol,
     $     alphaS,gammaS, deltaS,betaS,strainRateTension0,fc,fc0,
     $     strainRateCompression0,rate,ratioT,ratioC,rateFactorTension,
     $     rateFactorCompression,alpha,rateFactor,
     $     strainRateRatioCompression
      integer k
c     alpha                     --------------   alpha is the variable used in CDPM2U to evaluate contribution of compressive stresses to principal stress tensor  <Input>
c     totstrain(6)             --------------   array containing strain increments {xx,yy,zz,xy,yz,xz} <Input>
c     deltaTime                 --------------    time step  length <Input>
c     oldStrainRate                 --------------    strain rate  <Input/Output>
c     rateFactor                --------------   rateFactor used to incorporate impact effects on the constitutive law <Input>
c     fc                        --------------   concrete compressive strength <Input>
c     fc0                       --------------   concrete reference compressive strength as given in chapter 2.1.5.2 of Model Code 1990  <Input>
c     tol                       --------------   tolerance used to identify whether it is a tensile or compressive strain state <Input>
c     princDir(3,3)             --------------   matrix containing the principal directions of the strain rate tensor. Eigenvectors stored columnwise
c     princStrainRate(6)        --------------   array containing the 3 eigenvalues of the strain rate tensor. {rate1,rate2,rate3}
c     max                       --------------   max(tensile) principal strain rate
c     min                       --------------   min (compressive) principal strain rate
c     alphaS,betaS,gammaS,deltaS--------------   parameters for the formulation of the influence of the strain rate on the constitutive law given in eqs. 2.1-44a,2.1-44b, 2.1-46, 2.1-48a of the Model Code 1990 respectively
c     strainRateTension0        --------------   sig0Tension (see MC90)
c     strainRateCompression0    --------------   sig0Compression (see MC90)
c     ratioT                    --------------   max strainrate/sig0Tension (see MC90)
c     ratioC                    --------------   max strainrate/sig0Compression (see MC90)
c     rateFactorTension         --------------   contribution of tensile strain rates to the rate Factor
c     rateFactorCompression     --------------   contribution of compressive strain rates to the rate Factor
c     strainRate                 --------------    strain rate  

c     It is important that fc and fc0 are given in the same units
c     Therefore, we do not have a  default value for fc0 to avoid 
c     that users use for fc Pa but for the default value of fc0 MPa
      alphaS = 1. / ( 5. + 9. * fc / fc0 )
      gammaS = exp( ( 6.156 * alphaS - 2. ) * log(10.0) )
      deltaS = 1. / ( 1. + 8. * fc / fc0 )
      betaS = exp( ( 6. * deltaS - 2. ) * log(10.) )
      strainRateTension0 = 1.e-6
      strainRateCompression0 = -30.e-6
      rateFactorTension=1.
      rateFactorCompression=1.
      
      call cdpm2u_computePrincValues(strainRate,princStrainRate,1,
     $     princDir)
      
      max= -1.e-20
      min = 1.e20

      do k=1,3
         if (max .lt. princStrainRate(k)) then
            max=princStrainRate(k)           
         end if
         if (min .gt. princStrainRate(k)) then
            min=princStrainRate(k)
         end if
      enddo
      
      if ( 1. - alpha .gt.tol ) then 
         rate = max
      else
         rate =  min
      endif

      ratioT= rate/strainRateTension0
      if ( rate .lt. 30.e-6 ) then
         rateFactorTension = 1.
      else if ( 30.e-6 .lt. rate .and. rate .lt. 1.) then
         rateFactorTension = ratioT**deltaS
      else 
        rateFactorTension =  betaS * (ratioT **(1./3.))
      endif

      ratioC= rate/strainRateCompression0
      if ( rate .gt. -30.e-6 ) then
         rateFactorCompression = 1.
      else if (-30.e-6 .gt. rate .and. rate .gt. -30) then
         rateFactorCompression = ratioC**(1.026 * alphaS)        
      else 
         rateFactorCompression =  gammaS*(ratioC**(1./3.))
      endif

      rateFactor = ( 1. - alpha ) * rateFactorTension + 
     $     alpha * rateFactorCompression
      cdpm2u_computeRateFactor=rateFactor
      return
      end

      subroutine cdpm2u_computeDamage(omegaC,omegaT,strainRate,
     $     rateFactor,alpha,epsilonT,epsilonC,kappaDT,kappaDT1,kappaDT2,
     $     kappaDC,kappaDC1,kappaDC2,stress,deltaPlasticStrainNormT,
     $     tempKappaP,len,stressOld,oldAlpha,epsilon,epsilonNew)
c     Subroutine to perform the damage return. Both compressive and tensile damage variables are calculated.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag

      real*8 omegaC,omegaT,tempRateFactor,rateFactor,epsilonT,
     $     epsilonC,kappaDT,kappaDT1,kappaDT2,kappaDC,kappaDC1,kappaDC2,
     $     stress(6),deltaPlasticStrainNormT,deltaPlasticStrainNormC,
     $     tempKappaP,len,omegaOldC,omegaOldT,alpha,sigElastic,
     $     rhoElastic, thetaElastic,pHelp,help,tempEquivStrain,qhelp,
     $     tempEquivStrainT,tempEquivStrainC,fsT,fsC,ducMeas,Rs,
     $     strainRate(6),yieldTolDamage,cdpm2u_computeRateFactor,
     $     deltaPlasticStrainNormTNew,
     $     cdpm2u_computeDeltaPlasticStrainNormT,
     $     deltaPlasticStrainNormCNew,
     $     cdpm2u_computeDeltaPlasticStrainNormC,
     $     alphaZero,e0,rFunction,
     $     cdpm2u_computeDamageT,cdpm2u_computeDamageC,
     $     cdpm2u_computeEquivalentStrain,stressOld(6),
     $     minEquivStrain,oldAlpha,epsilon,epsilonNew
      integer unloadingFlag,step1flag
c     omegaC          --------------- compressive damage variable                        <Input/Output>
c     omegaT          --------------- tensile damage variable                            <Input/Output>
c     totstrain(6)    --------------- array containing the rate of the strain tensor     <Input>
c     strainRate      ---------------  strain rate     <Input/Output>
c     deltaTime       ---------------  time step length     <Input>
c     rateFactor      --------------- variable to incorporate impact effects on the constitutive law  <Input/Output>
c     alpha           --------------- variable showing the contribution of damage on the strains      <Input/Output>
c     epsilonT        --------------- tensile equivalent strain                <Input/Output>
c     epsilonC        --------------- compressive equivalent strain            <Input/Output>
c     epsilon         ---------------  old equivalent strain (no rate influnce)            <Input/Output>
c     kappaDT         --------------- history parameter kappaDT                <Input/Output>
c     kappaDT1        --------------- history parameter kappaDT1               <Input/Output>
c     kappaDT2        --------------- history parameter kappaDT2               <Input/Output>
c     kappaDC         --------------- history parameter kappaDC                <Input/Output>
c     kappaDC1        --------------- history parameter kappaDC1               <Input/Output>
c     kappaDC2        --------------- history parameter kappaDC2               <Input/Output>
c     stress(6)       --------------- stress tensor                            <Input>
c     deltaPlasticStrainNormT ------- norm of the increment of plastic strains <Input>
c     tempKappaP      --------------- cummulative plastic strain (kappaP)      <Input>
c     deltaElStrain(6)  ---------------  increment of the elastic strain vector  <Input>
c     oldAlpha        ---------------  old alpha  <Input>
c     len             --------------- characteristic length used to combine damage law with the crack-band approach <Input>
c     deltaPlasticStrainNormTNew ---- norm of the increment of the cummulative plastic strain used in the damage algorithm for tension
c     cdpm2u_computeDeltaPlasticStrainNormT- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the tensile damage variable
c     deltaPlasticStrainNormCNew ---- norm of the increment of the cummulative plastic strains used in the damage algorithm for compression
c     cdpm2u_computeDeltaPlasticStrainNormC- function used to calculate the increment of the cummulative plastic strain kappa used for the calculation of the compressive damage variable
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     alphaZero       --------------- parameter used for the calculation of the damage ductility measure based on eqs.(56-57) of the IJSS paper by Grassl et al.
c     e0              --------------- parameter equal to ft/E
c     cdpm2u_computeDamageT,cdpm2u_computeDamageC - functions used to calculate tensile and compressive damage variables respectively
c     yieldTolDamage  --------------- tolerance used in the damage algorithm for cased where Hp=0
c     Rs              --------------- variable used in the calculation of the ductility measure for damage given in eq.(57) of the IJSS paper by Grassl et al.
c     xs              --------------- ductility measure for damage given in eq.(56) of the IJSS paper by Grassl et al.
c     fsT,fsC         --------------- loading function for the tensile and compressive damage variables
c     cdpm2u_computeRateFactor ------------- function used for the calculation of the rate factor to include impact effects on the constitutive law
c     omegaOldT,omegaOldC ----------- old (previous step) tensile and compressive damage variables
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression)
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     minEquivStrain    -------------- when unloading is occuring corresponds to the minEquivStrain before reloading occurs
c     step1flag    -------------- flag indicating whether we are at the first step
c                                        =1 it is analysis 1st step
c                                        =0  it not analysis 1st step

      e0=ft/ym
      yieldTolDamage=gTol*10.0
      omegaOldC=omegaC
      omegaOldT=omegaT
      deltaPlasticStrainNormC=deltaPlasticStrainNormT
      if (rateFactor .eq. 0.) then
         step1flag=1
         rateFactor=1.
      else
         step1flag=0
      end if
      call cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,tempEquivStrain,epsilon,ym,pr,
     $     gtol)

c-------------------Compute tensile and compressive equivalent strains-------------------------------------
      if (strrateflg .eq. 1.0 .and. omegaC .eq. 0. .and. 
     $     omegaT.eq. 0. ) then
         tempRateFactor=cdpm2u_computeRateFactor(alpha,strainrate,
     $        fc,fc0,gTol)
      else 
         tempRateFactor=rateFactor
      end if
	  
      tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/
     $     tempRateFactor
      
      if (unloadingFlag .eq. 0) then
         tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $        tempRateFactor
      else
          tempEquivStrainC=epsilonC+ oldAlpha*(minEquivStrain-epsilon)/
     $        tempRateFactor + alpha*(tempEquivStrain-minEquivStrain)/
     $         tempRateFactor
      end if
   
c     Note rate factor is calculated only once at the onset of damage
      if ( ( tempEquivStrainT .gt. e0 .or. tempEquivStrainC .gt. e0
     $     ) .and. ( ( omegaT .eq. 0. ) .and. (omegaC .eq. 0. ) ).and.
     $     strrateflg .eq. 1.0 .and. step1flag .ne. 1) then
         tempEquivStrainT=epsilonT+(tempEquivStrain-epsilon)/rateFactor
         if (unloadingFlag .eq. 0) then
            tempEquivStrainC=epsilonC+(tempEquivStrain-epsilon)*alpha/
     $           rateFactor
         else
            tempEquivStrainC=epsilonC+oldAlpha*(minEquivStrain-
     $           epsilon)/rateFactor+(tempEquivStrain-minEquivStrain)*
     $           alpha/rateFactor 
         end if
      else
         rateFactor = tempRateFactor 
      endif
      
      fsT = (tempEquivStrainT - kappaDT)/e0
      fsC = (tempEquivStrainC - kappaDC)/e0

      epsilonNew=tempEquivStrain
      epsilonT=tempEquivStrainT
      epsilonC=tempEquivStrainC
c -------------------- Compute Ductility Measure Damage --------------------------------------------------
      call cdpm2u_computeTrialCoordinates(stress, sigElastic,
     $rhoElastic,thetaElastic)
      Rs = 0.
      alphaZero= 1./sqrt(6.) 
      if ( sigElastic .lt. 0. ) then
         if ( rhoElastic .gt. 1.e-16 ) then
            Rs = -sigElastic /(alphaZero*rhoElastic)
         else 
            Rs = -sigElastic * 1.e16 / alphaZero
         endif
      else 
         Rs = 0.
      endif
      ducMeas = 1. + ( as - 1. ) * Rs
      
c----- Check which damage surfaces (tensile/compressive) are active --------------------------------------
      if (fsT .lt. -yieldTolDamage .and. 
     $     fsC .lt. -yieldTolDamage) then
c     no increase of the damage variables required
      else if (  fsT .ge. -yieldTolDamage .and. 
     $        fsC .lt. -yieldTolDamage) then
c     only tensile damage surface active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT, deltaPlasticStrainNormT,kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / ducMeas / rateFactor
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / 
     $        ducMeas
         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT)
      else if (  fsT .lt.  -yieldTolDamage .and. 
     $        fsC .ge.  -yieldTolDamage) then
c     only compressive damage surface active
         deltaPlasticStrainNormCNew =
     $        cdpm2u_computeDeltaPlasticStrainNormC(alpha, 
     $        tempEquivStrainC,
     $        deltaPlasticStrainNormC, kappaDC,rhoElastic,
     $        tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas*rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC

         omegaC = 
     $        cdpm2u_computeDamageC(kappaDC, 
     $        kappaDC1, kappaDC2, omegaOldC)
      else if (  fsT .ge.-yieldTolDamage  .and. 
     $        fsC .ge.-yieldTolDamage ) then
c Both compressive and tensile damage surfaces are active
         deltaPlasticStrainNormTNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormT(
     $        tempEquivStrainT,deltaPlasticStrainNormT, kappaDT,
     $        e0,yieldTolDamage)
         kappaDT1 = kappaDT1 + 
     $        deltaPlasticStrainNormTNew / (ducMeas * rateFactor)
         kappaDT2 = kappaDT2 + ( tempEquivStrainT - kappaDT) / ducMeas

         
         kappaDT= tempEquivStrainT

         omegaT = cdpm2u_computeDamageT(kappaDT, kappaDT1, kappaDT2, 
     $        len, omegaOldT)

         deltaPlasticStrainNormCNew = 
     $        cdpm2u_computeDeltaPlasticStrainNormC(
     $        alpha,tempEquivStrainC,deltaPlasticStrainNormC, 
     $        kappaDC,rhoElastic,tempKappaP,yieldTolDamage)
         kappaDC1 = kappaDC1 + 
     $        deltaPlasticStrainNormCNew /( ducMeas * rateFactor)
         kappaDC2 = kappaDC2 + ( tempEquivStrainC - kappaDC) / 
     $        ducMeas
         
         kappaDC= tempEquivStrainC
         omegaC = cdpm2u_computeDamageC(kappaDC, kappaDC1, 
     $        kappaDC2, omegaOldC)
      endif
    
      return
      end

      subroutine cdpm2u_checkForUnAndReloading(stress,stressOld,
     $     unloadingFlag,minEquivStrain,equivStrainNew,equivStrainOld,
     $     ym,pr,gtol)
c     Function to check if there is unloading and reloading occuring within one step (usually happens during cyclic loading). If such a process is happening the algorithm returns a flag and an approximate value of the minimum equivalent strain. Moreover it returns the current equivalent strain.

      real*8 stress(6),minEquivStrain,stress1(6),ym,pr,
     $     deltaStress(6),stressPlus(6),stressMinus(6),
     $     equivStrainOld,equivStrainNew,equivStrain1,stressOld(6),
     $     equivStrainPlus,equivStrainMinus,sigV,rho,theta,
     $     cdpm2u_computeEquivalentStrain,gtol
      integer i,j,unloadingFlag
c     stress(6)   ------------- effective stress tensor(no damage included) <Input>
c     stressOld(6) ---------- effective stress vector in  previous step <Input>
c     equivStrainOld ---------- equivalent strain of the previous loading step <input>
c     ym             ---------- Young's modulus <Input>
c     pr             ---------- Poisson's ratio <Input>
c     unloadingFlag    -------------- flag indicating whether unloading and reloading is occuring during the current step (e.g. transition from tension to compression) <Output>
c                                        =1 unloading and reloading occur within a step
c                                        =0  no unloading and reloading occur within a step
c     equivStrainNew ---------- equivalent strain of the current loading step  <Output>
c     minEquivStrain   -------- minimum equivalent strain <Output>
c     stressPlus    ----------   stress vector representing previous step's elastic strain vector plus 0.99*deltaStrain
c     stressMinus   ----------   stress vector representing current step's elastic strain vector minus 0.01*deltaStrain
c     sigV          ---------- volumentric stress
c     rho           ---------- deviatoric stress
c     theta         ---------- Lode angle
c     equivStrain1  ---------- equivalent strain of various strain vectors
c     computeEquivalentStrain -------- function to calculate equivalent strain according to eq. 37 if the IJSS paper by Grassl et al

      call cdpm2u_computeTrialCoordinates(stress, sigV, rho,theta)
      equivStrainNew= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
    
      do i=1,6
        deltaStress(i)=stress(i)-stressOld(i)
      enddo
      do i=1,6
         stressPlus(i)=stressOld(i)+0.01*deltaStress(i)
         stressMinus(i)=stressOld(i)+0.99*deltaStress(i)
      enddo

      call cdpm2u_computeTrialCoordinates(stressPlus, sigV, rho,theta)
      equivStrainPlus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      call cdpm2u_computeTrialCoordinates(stressMinus, sigV, rho,theta)
      equivStrainMinus= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 

      unloadingFlag=0
      minEquivStrain=equivStrainOld
      if ( (equivStrainPlus .lt. equivStrainOld .and.
     $     equivStrainMinus .lt. equivStrainNew) .and. 
     $     (abs(equivStrainPlus - equivStrainOld) .gt. gtol/10. .and. 
     $     abs(equivStrainMinus - equivStrainNew) .gt. gtol/10.)) then
         unloadingFlag=1
         write(*,*) '*** Unloading and reloading occurs.'
         write(*,*) 'within a single step. Subincrementation performed' 
         do i=1,100
            do j=1,6
               stress1(j)=stressOld(j)+deltaStress(j)*FLOAT(i)/100.0 
            enddo
            call cdpm2u_computeTrialCoordinates(stress1,sigV,rho,theta)
            equivStrain1= cdpm2u_computeEquivalentStrain(sigV,rho,theta) 
            if( equivStrain1 .le. minEquivStrain) then
               minEquivStrain=equivStrain1
            else
               goto 625
            end if
         enddo
      end if
      
 625  continue
      return
      end
      
      real*8 function  cdpm2u_computeEquivalentStrain(sigElastic,
     $     rhoElastic, thetaElastic)
c     Function to calculate equivalent strain used in the damage part according to eq. 37 of the IJSS paper by Grassl et al.

      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8  answer, rFunction,thetaElastic,rhoElastic,
     $     sigElastic, pHelp,qHelp,help,e0
c     sigElastic, rhoElastic, thetaElastic ---- volumetric and deviatoric stresses and Lode angle of the effective stress tensor 
c     pHelp,help,qhelp -------------- variables used in the calculation of the equivalent strain
c     rFunction       --------------- parameter given by eq.(19) of the IJSS paper by Grassl et al.      
c     e0              --------------- parameter equal to ft/E
      e0=ft/ym
      rFunction = ( 4. * ( 1. - ecc**2. ) * cos(thetaElastic)**2. +
     $     ( 2. * ecc - 1. )**2.  ) /
     $     ( 2. * ( 1. - ecc**2. ) * cos(thetaElastic) +
     $     ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc**2. ) *
     $     cos(thetaElastic)**2.+ 5. * ecc**2. - 4. * ecc) )
      pHelp = -m0 * ( rhoElastic * rFunction / ( sqrt(6.) * fc ) + 
     $     sigElastic / fc )
      qHelp = -3. / 2. * (rhoElastic** 2.) / (fc**2.)
      help = -0.5 * pHelp + sqrt((pHelp** 2.) / 4. - qHelp)
c     negative help values are not of interest and create problems since we compute the square root 
      answer = 0.
      if ( help .gt. 0. ) then
         answer= help * e0
      endif
      cdpm2u_computeEquivalentStrain=answer
      return 
      end
      
      real*8 function  cdpm2u_computeDamageT(kappa,kappaOne,kappaTwo,le, 
     $      omegaOld)
c     Function to calculate damage due to tension.
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 kappa, kappaOne, kappaTwo, omegaOld,residual,le,
     $     residualDerivative,omega,tol,help,e0,yieldTolDamage
      integer iter,newtonIter
c     kappa           ----------------------- damage history parameter kappaDT <input>
c     kappaOne        ----------------------- damage history parameter kappaDT1 <input>
c     kappaTwo        ----------------------- damage history parameter kappaDT2 <input
c     omegaOld        ----------------------- old damage variable (previous step) omegaT <input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     residualDerivative--------------------- residual used to calculate with N-R iterative procedure the omegaT for exponential softening law
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaT
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaT for exponential softening law
c     help            -----------------------  help variable used in the billinear softening law
c     iter            ----------------------- counter of performed N-R iterations in the exponential softening law
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaT in the exponential softening law
c     yieldTolDamage  ----------------------- tolerance used in the damage algorithm if Hp=0


      newtonIter=100
      tol=gTol/100.0
      e0=ft/ym
      yieldTolDamage=gTol*10.0
      if ( kappa  .gt. e0*(1-yieldTolDamage) ) then
        if ( type .eq. 0. ) then
c     Linear damage law
           omega = ( ym * kappa * wf - ft * wf + ft * kappaOne * le ) /
     $          ( ym * kappa * wf - ft * le * kappaTwo )
        else if ( type .eq. 1. ) then
c     Billinear damage law
            omega = ( ym * kappa * wf1 - ft * wf1 - ( ft1 - ft ) * 
     $          kappaOne * le ) /( ym * kappa * wf1 + 
     $          ( ft1 - ft ) * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo
            if ( help .ge. 0. .and. help .lt. wf1 ) then
               goto 185
            endif

            omega = ( ym * kappa * ( wf - wf1 ) - ft1 * ( wf - wf1 ) +
     $           ft1 * kappaOne * le  - ft1 * wf1 ) / ( ym * kappa * 
     $           ( wf - wf1 )  - ft1 * le * kappaTwo )
            help = le * kappaOne + le * omega * kappaTwo

            if ( help .gt. wf1 .and. help .lt. wf ) then
               goto 185
            endif
         else if ( type .eq. 2. ) then
c     Exponential: Iterative solution with N-R procedure
            omega = 1.
            residual=0.
            residualDerivative = 0.
            iter = 0

 135        continue
            iter=iter+1
            residual = ( 1 - omega ) * ym * kappa - ft *
     $           exp(-le * ( omega * kappaTwo + kappaOne ) / wf)
            residualDerivative = -ym * kappa + ft * le * 
     $           kappaTwo / wf * exp(-le * ( omega * kappaTwo + 
     $           kappaOne ) / wf)
            omega =omega- residual / residualDerivative;
            if ( iter .gt. newtonIter ) then
               write(*,*) '*** Algorithm for tensile damage-
     $No convergence reached after 100 iterations ***'
               stop
            end if
            if (abs(residual/ft) .ge. 1.e-8) then
               goto 135
            end if
         end if
      else 
         omega = 0.;
      endif
      
 185  continue
 
      if ( omega .gt. 1. ) then
         omega = 1.
      endif
       
      if ( omega .lt. 0. .or. omega .lt. omegaOld) then
         omega=omegaOld
      endif
      cdpm2u_computeDamageT= omega
      return 
      end  

      real*8 function cdpm2u_computeDamageC(kappa, kappaOne, 
     $     kappaTwo, omegaOld)
c     Function to calculate damage due to compression
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      
      real*8 kappa, kappaOne, kappaTwo, omegaOld,residual,dResidualDOmega,
     $     exponent,omega,tol,kappaDC,e0,yieldTolDamage,errorOld
      integer nite,newtonIter
c     kappa           ----------------------- damage history parameter kappaDC <Input>
c     kappaOne        ----------------------- damage history parameter kappaDC1 <Input>
c     kappaTwo        ----------------------- damage history parameter kappaDC2 <Input
c     omegaOld        ----------------------- old damage variable (previous step) omegaC <Input>
c     residual        ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     dResidualDOmega ----------------------- residual used to calculate with N-R iterative procedure the omegaC
c     exponent        ----------------------- exponent used in the formulation of the damage law (in the current version assumed =1.)
c     omega           ----------------------- new damage variable calculated based on the input of the current step omegaC
c     tol             ----------------------- tolerance used in N-R procedure for the calculation of omegaC
c     nite            ----------------------- counter of performed N-R iterations
c     newtonIter      ----------------------- value of the max allowed N-R iterations for the calculation of omegaC
      if (isoflag .eq. 1.0) then
         omega=0.0
         kappaOne=0.0
         kappaTwo=0.0
         kappa=0.0
         goto 188
      end if

c      newtonIter=100
      newtonIter=200
c      tol=gTol/100.0
      tol=gTol
      yieldTolDamage=gTol*10.

      omega=1.
      nite = 0
      residual = 0.
      dResidualDOmega = 0.
      exponent = 1.
      e0=ft/ym
      if ( kappa .gt. e0*(1-yieldTolDamage)) then
 187     continue 
         nite=nite+1;
         residual =  ( 1. - omega ) * ym * kappa - ft * exp( - ( 
     $        kappaOne + omega * kappaTwo ) / efc )
         dResidualDOmega =-ym * kappa + ft * kappaTwo / efc * exp( -( 
     $        kappaOne + omega * kappaTwo ) / efc )
         omega = omega- residual / dResidualDOmega
         errorOld = residual/ft
         if (omega .lt. 0.) then
            omega = 0.
            goto 188
         end if 
         if ( nite .gt. newtonIter ) then
            if (residual .lt. 0.) then
               omega = omegaOld
               goto 188
            end if
            write(*,*) '*** Algorithm for compressive damage-
     $No convergence reached.'
c Set omega=omegaOld =', omegaOld,' ***'
c            omega = omegaOld
            write(*,*) '   (current) residual/ft = ',residual/ft
            write(*,*) '   (prev.)   residual/ft = ',errorOld
            write(*,*) '   omega                 = ',omega
            write(*,*) '   omegaOld              = ',omegaOld
            stop
         endif
         if( abs(residual/ft) .ge. tol )   goto 187
      else 
         omega = 0.
      endif
	   
 188  continue

      if ( omega .gt. 1 ) then
         omega = 1.
      endif
       
      if ( omega .lt. 0. .or. omega .lt. omegaOld ) then
         omega = omegaOld
      endif
       
      cdpm2u_computeDamageC=omega
      return
      end

      


      real*8 function  cdpm2u_computeDeltaPlasticStrainNormT(tempKappaD,
     $     plastStrNorm, kappaD,e0,yieldTolDamage)      
c     Function returning the norm of the increment of the plastic strain tensor.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch

      real*8 e0,yieldTolDamage,answer,tempKappaD,plastStrNorm, kappaD,
     $     factor
c     tempKappaD      -------------------    temporary (current) KappaDt                   <Input>
c     kappaD          -------------------    old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------    plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------    variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------    tolerance used when Hp=0                      <Input>
c     factor          -------------------    factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     answer          -------------------    calculated norm             
      factor = 0.
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0  * ( 1. - yieldTolDamage )) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer = plastStrNorm*factor
      else 
         answer=plastStrNorm 
      end if
      cdpm2u_computeDeltaPlasticStrainNormT=answer
      return 
      end

      real*8 function  cdpm2u_computeDeltaPlasticStrainNormC( 
     $     alpha,tempKappaD,
     $     plastStrNorm, kappaD,rho,tempKappa,yieldTolDamage)
c     Function returning the norm of the increment of the plastic strain tensor multiplied by alphaC and betaC according to eq. (48) of the IJSS paper by Grassl et al.
c     Special treatment is applied during transition from hardening (pre-peak) to the post-peak branch      
      real*8 ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      common/cdpmc/ym,pr,ecc,qh0,ft,fc,hp,ah,bh,ch,dh,as,df,fc0,type,bs,
     1     wf,wf1,efc,ft1,strrateflg,failflg,m0,gTol,isoflag,printflag
      real*8   tempKappaD,plastStrNorm, kappaD,factor,qh2,alpha,rho,e0,
     $     cdpm2u_qh2fun,extraFactor,tempKappa,kappa,
     $     yieldTolDamage,answer

c     tempKappaD      -------------------  temporary (current) KappaDt                   <Input>
c     kappaD          -------------------  old (previous step) KappaDt                   <Input>
c     plastStrNorm    -------------------  plastic strain incremement norm (epNew-epOld) <Input>
c     e0              -------------------  variable equal to ft/E                        <Input>
c     yieldTolDamage  -------------------  tolerance used when Hp=0                      <Input>
c     rho             -------------------  deviatoric stress                             <Input>
c     alpha           -------------------  variable calculated based on eq. 46 of the IJSS paper by P. Grassl et al.    <Input>
c     tempKappa       -------------------  temporary (current step) cummulative plastic strain (kappaP)   <Input>
c     factor          -------------------  factor to multiply the calculated norm during transition from the hardening(pre-peak) to softening (post-peak)   
c     qh2             -------------------  variables containing the results of the hardening function 
c     cdpm2u_qh2fun          -------------------  functions to calculate the hardening functions  given in eq. (31) in IJSS paper by Grassl et al.
c     extraFactor     -------------------  variable betaC given in eq. (50) of the IJSS paper by Grassl et al.
c     answer          -------------------  calculated answer         
      e0=ft/ym
      if ( tempKappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         answer = 0.
      else if ( tempKappaD .gt. e0 * ( 1. - yieldTolDamage ) .and.
     $        kappaD .lt. e0 * ( 1. - yieldTolDamage ) ) then
         factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) )
         answer=plastStrNorm*factor
      else 
         answer=plastStrNorm
      end if
      
      qh2=cdpm2u_qh2fun(tempKappa,hp)
      
      if (rho<1.e-16) then
         extraFactor =ft * qh2 * sqrt(2. / 3.) / 1.e-16 / sqrt( 1. + 
     $        2.*  (df** 2.) )
      else 
         extraFactor =ft * qh2 * sqrt(2. / 3.) / rho / sqrt( 1. + 
     $        2.* (df** 2.) )
      endif
      answer=answer*extraFactor*alpha
      cdpm2u_computeDeltaPlasticStrainNormC=answer
      return
      end
      

c ---------------------------------------------General functions--------------------------------

      subroutine cdpm2u_computePrincValues(ax,w,flag,z)
c     Subroutine to solve the eigenvalues and eigenvectors of real*8 symmetric matrix by jacobi method.
c     Eigenvalues are stored in descending order (from max to min) in array w and eigenvectors
c     are stored columnwise in matrix z.

      real*8 ax,w,z,matrix
      integer flag
      dimension ax(6),w(3), z(3,3),matrix(3,3)
      real*8 ssum, aa, co, si, tt, tol, sum, aij, aji;
      integer ite, i, j, k, ih
      
      real*8 swap
      integer ii,jj,kk
c     ax(6)    ---------------------            Input tensor. It is a 6x1 array and components are given 
c                                               in order xx,yy,zz,xy,yz,xz  <Input>
c     flag     ---------------------            Flag denoting whether it is a stress or a strain tensor in 
c                                                order to convert voigt shear strains(gammas) to tensorial
c                                                shear strains(epsilons).
c                                                = 0: Stress tensor
c                                                = 1: Strain tensor
c     w(3)     --------------------             3x1 Array containing eigenvectors of ax stored in 
c                                               descending order(max to min)  <Output>
c     z(3,3)   --------------------             Matrix containing all eigenvectors stored columnwise <Output>
c     matrix(3,3) -----------------             Stress/Strain tensor that is being processed
c     tol      --------------------             Tolerance of the algorithm (default 10 figures,1e-10)
c     ite,i,j,k,ih ----------------             Counters used in various loops
c     sum,ssum --------------------             Summation variables used in various addition procedures
c     aa,co,si,tt,aij,aji ---------             Variables used in the algorithm 

      tol=1.e-10
c     Reconstruct tensor based on given array
      matrix(1,1)=ax(1)
      matrix(2,2)=ax(2)
      matrix(3,3)=ax(3)
      if (flag .eq. 0) then
         matrix(2,1)=ax(4)
         matrix(1,2)=ax(4)
         matrix(3,1)=ax(6)
         matrix(1,3)=ax(6)
         matrix(3,2)=ax(5)
         matrix(2,3)=ax(5)
      else
         matrix(2,1)=ax(4)/2.
         matrix(1,2)=ax(4)/2.
         matrix(3,1)=ax(6)/2.
         matrix(1,3)=ax(6)/2.
         matrix(3,2)=ax(5)/2.
         matrix(2,3)=ax(5)/2.
      endif
      
c     Initialise w,z and check if zero stress state
      do i=1,3
         w(i) = matrix(i, i)
      enddo
      sum=0.
      do i=1,3
         do j=1,3
            sum =sum+ abs( matrix(i, j) )
            z(i, j) = 0.0
         enddo
         z(i, i) = 1.0
      enddo
      if ( sum .le. 0.0 ) then
         goto 900
      endif
            
c     Reduce to matrix diagonal
      ite=0

 272  continue
      ssum = 0.0
      do j=2,3
         ih = j - 1
         do i=1,ih
            if ( abs( matrix(i, j) ) / sum  .gt. tol ) then
               ssum =ssum+ abs( matrix(i, j) )
c     CALCULATE ROTATION ANGLE
               aa = atan2( matrix(i, j) * 2.0, w(i) - w(j) ) /  2.0
               si = sin(aa)
               co = cos(aa)
               
c     MODIFY "I" AND "J" COLUMNS OF "matrix" and "z"
               do k=1,i-1
                  tt = matrix(k, i)
                  matrix(k, i) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (i,i)
               tt = w(i)
               w(i) = co * tt + si *matrix(i, j)
               aij = -si * tt + co *matrix(i, j)
               tt = z(i, i)
               z(i, i) = co * tt + si *z(i, j)
               z(i, j) = -si * tt + co *z(i, j)
               
               do k=i+1,j-1
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(k, j)
                  matrix(k, j) = -si * tt + co *matrix(k, j)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     diagonal term (j,j)
               tt = matrix(i, j)
               aji = co * tt + si *w(j)
               w(j) = -si * tt + co *w(j)
               
               tt = z(j, i)
               z(j, i) = co * tt + si *z(j, j)
               z(j, j) = -si * tt + co *z(j, j)
               
               do k=j+1,3
                  tt = matrix(i, k)
                  matrix(i, k) = co * tt + si *matrix(j, k)
                  matrix(j, k) = -si * tt + co *matrix(j, k)
                  tt = z(k, i)
                  z(k, i) = co * tt + si *z(k, j)
                  z(k, j) = -si * tt + co *z(k, j)
               enddo
c     MODIFY DIAGONAL TERMS
               w(i) = co * w(i) + si * aji
               w(j) = -si * aij + co *w(j)
               matrix(i, j) = 0.0
            else 
c     matrix(I,J) MADE ZERO BY ROTATION
            endif
         enddo
      enddo
      
      ite=ite+1
c     CHECK FOR CONVERGENCE
      if ( ite .gt. 50 ) then
         write(*,*)  '*** Compute principal values.
     $Too many iterations! ***'
         stop
      endif
      
      if ( abs(ssum) / sum .gt. tol  ) goto 272
      
      do ii=1,2
         do jj=1,2
            if ( w(jj + 1) > w(jj) ) then
c     swap eigenvalues and eigenvectors
               swap = w(jj + 1);
               w(jj + 1) = w(jj);
               w(jj) = swap;
               do kk=1,3
                  swap = z(kk, jj + 1);
                  z(kk, jj + 1) = z(kk, jj);
                  z(kk, jj) = swap;
               enddo
            endif
         enddo
      enddo
 900  continue
      return
      end
      
      subroutine cdpm2u_computeInverseJac(jacinv,jac,error)
c     Subroutine to calculate the inverse of a 4x4 matrix
      real*8 jacinv(4,4),jac(4,4),tmp(4,4)
      real*8  piv, linkomb,dtol
      integer i,j,k,error
c     jac(4,4)    ---------  matrix whose inverse will be calculated <Input>
c     jacinv(4,4) ---------  inverse matrix  <Output>
c     error       --------- integer showing whether solution has converged <Output>
c                            = 0 solution has converged
c                            =-1 solution has not converged
c     tmp(4,4)    --------- temporary matrix
c     piv,lincomb --------- variables used in gaussian elimination
c     i,j,k       --------- counters
c     dtol        --------- tolerance of the algorithm

      dtol=1.e-20           
      do i=1,4
         do j=1,4
            tmp(i,j)=jac(i,j)
            jacinv(i,j)=0.
         enddo
         jacinv(i,i)=1.
      enddo
      
      do i=1,3
         piv = tmp(i, i)
         if (abs(piv) .lt. dtol) then
            error=-1
            goto 452
         endif

         do j=i+1,4
            linkomb = tmp(j, i) / tmp(i, i)
            do k=i,4
               tmp(j, k) = tmp(j, k) - tmp(i, k) * linkomb
            enddo
            do k=1,4
               jacinv(j, k) =jacinv(j, k)-jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=4,2,-1
         piv = tmp(i, i)
         do j=i-1,1,-1
            linkomb = tmp(j, i) / piv
            do k=i,1,-1
               tmp(j, k) =  tmp(j, k) -tmp(i, k) * linkomb
            enddo
            do k=4,1,-1
               jacinv(j, k) =jacinv(j, k) - jacinv(i, k) * linkomb
            enddo
         enddo
      enddo
      
      do i=1,4
         do j=1,4
            jacinv(i, j) = jacinv(i, j) / tmp(i, i)
         enddo
      enddo 
      error=0
 452  continue
      return
      end
      
      
      
      subroutine cdpm2u_transformStressVectorTo(stress,princDir,
     $     princStress)
c     Rotates principal stress tensor back to the original Coordinate system. It calculates transformation 
c     matrix  and then multiplies with principal stress tensor.

      real*8 stress(6),princDir(3,3),princStress(3),transpose(3,3),
     $     transformationMtrx(6,6),sum,princStressTensor(6)
      integer i,j
c     stress (6)             ----------------- stress vector <Output>
c     princDir(3,3)          ----------------- matrix containing eigenvectors stored columnwise <Input>
c     princStress(3)         ----------------- stress vector containing principal stresses =
c                                              {sigma1,sigma2,sigma3}  <Input>
c     princStressTensor(6)   ----------------- stress vector at principal axis coordinate system 
c                                              ={sigma1,sigma2,sigma3,0,0,0}
c     transpose(3,3)         ----------------- matrix containing eigenvectors stored rowise
c     sum                    ----------------- variable used in summation during matrix multiplication
c     transformationMtrx(6,6)----------------- transformation matrix used to transform principal stress
c                                              vector to stress vector in original CS
c     i and j                -----------------  integers used as counters


      do i=1,3
         princStressTensor(i)= princStress(i)
         princStressTensor(i+3)= 0.
         do j=1,3
            transpose(i,j)=princDir(j,i)
         enddo
      enddo
      
      
      transformationMtrx(1,1)=transpose(1,1)*transpose(1,1)
      transformationMtrx(1,2)=transpose(2,1)*transpose(2,1)
      transformationMtrx(1,3)=transpose(3,1)*transpose(3,1)
      transformationMtrx(1,4)=2.*transpose(1,1)*transpose(2,1)
      transformationMtrx(1,5)=2.*transpose(2,1)*transpose(3,1)
      transformationMtrx(1,6)=2.*transpose(1,1)*transpose(3,1)

      transformationMtrx(2,1)=transpose(1,2)*transpose(1,2)
      transformationMtrx(2,2)=transpose(2,2)*transpose(2,2)
      transformationMtrx(2,3)=transpose(3,2)*transpose(3,2)
      transformationMtrx(2,4)=2.*transpose(1,2)*transpose(2,2)
      transformationMtrx(2,5)=2.*transpose(2,2)*transpose(3,2)
      transformationMtrx(2,6)=2.*transpose(1,2)*transpose(3,2)

      transformationMtrx(3,1)=transpose(1,3)*transpose(1,3)
      transformationMtrx(3,2)=transpose(2,3)*transpose(2,3)
      transformationMtrx(3,3)=transpose(3,3)*transpose(3,3)
      transformationMtrx(3,4)=2.*transpose(1,3)*transpose(2,3)
      transformationMtrx(3,5)=2.*transpose(2,3)*transpose(3,3)
      transformationMtrx(3,6)=2.*transpose(1,3)*transpose(3,3)

      transformationMtrx(4,1)=transpose(1,1)*transpose(1,2)
      transformationMtrx(4,2)=transpose(2,1)*transpose(2,2)
      transformationMtrx(4,3)=transpose(3,1)*transpose(3,2)
      transformationMtrx(4,4)=transpose(1,1)*transpose(2,2)+
     $     transpose(2,1)*transpose(1,2)
      transformationMtrx(4,5)=transpose(2,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(2,2)
      transformationMtrx(4,6)=transpose(1,1)*transpose(3,2)+
     $     transpose(3,1)*transpose(1,2)

      transformationMtrx(5,1)=transpose(1,2)*transpose(1,3)
      transformationMtrx(5,2)=transpose(2,2)*transpose(2,3)
      transformationMtrx(5,3)=transpose(3,2)*transpose(3,3)
      transformationMtrx(5,4)=transpose(1,2)*transpose(2,3)+
     $     transpose(2,2)*transpose(1,3)
      transformationMtrx(5,5)=transpose(2,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(2,3)
      transformationMtrx(5,6)=transpose(1,2)*transpose(3,3)+
     $     transpose(3,2)*transpose(1,3)
    

      transformationMtrx(6,1)=transpose(1,1)*transpose(1,3)
      transformationMtrx(6,2)=transpose(2,1)*transpose(2,3)
      transformationMtrx(6,3)=transpose(3,1)*transpose(3,3)
      transformationMtrx(6,4)=transpose(1,1)*transpose(2,3)+
     $     transpose(2,1)*transpose(1,3)
      transformationMtrx(6,5)=transpose(2,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(2,3)
      transformationMtrx(6,6)=transpose(1,1)*transpose(3,3)+
     $     transpose(3,1)*transpose(1,3)

 
      do i=1,6
         sum=0.
         do j=1,6
            sum=sum+  transformationMtrx(i,j)*princStressTensor(j)
         enddo
         stress(i)=sum
      enddo
c
      return
      end

	  
      subroutine cdpm2u_computeMatrixDeterminant(dim,mat,det)
c     Subroutine to calculate the determinant of a square matrix
      real*8 dim,iq,xx1,xx2,yy1,yy2
      real*8 mat(dim,dim), tmp(dim,2*dim-1)
      integer i,j
c     dim              --------- dimension of a square matrix <Input>
c     mat(dim,dim)     --------- square matrix  <Output>
c     tmp(dim,2*dim-1) --------- temporary matrix
c     det              --------- determinant of matrix
c     i,j,k       --------- counters

      do i=1,dim
         do j=1,dim
            tmp(i,j)=mat(i,j)
         enddo
      enddo
	  
      do i=1,dim
         do j=1,dim-1
            tmp(i,j+dim)=mat(i,j)
         enddo
      enddo
	  
      iq=0
      xx1=1.0d0
      xx2=0.0d0
      do i=1,dim
         xx1=1.0d0
         do j=1,dim
            xx1=xx1*tmp(j,j+iq)
         enddo
      xx2=xx2+xx1
      iq=iq+1
      enddo
	  
      iq=0
      yy1=1.0d0
      yy2=0.0d0
      do i=1,dim
         yy1=1.0d0
         iq=i+dim-1
         do j=1,dim
            yy1=yy1*tmp(j,iq)
            iq=iq-1
         enddo
      yy2=yy2+yy1
      enddo
	  
      det=xx2-yy2
	  
      return
      end


	  
	  
      subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, 
     &  rot_num )

c*********************************************************************72
c
cc JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
c
c  Discussion:
c
c    This function computes the eigenvalues and eigenvectors of a
c    real symmetric matrix, using Rutishauser's modfications of the classical
c    Jacobi rotation method with threshold pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 September 2013
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the matrix, which must be square, real,
c    and symmetric.
c
c    Input, integer IT_MAX, the maximum number of iterations.
c
c    Output, double precision V(N,N), the matrix of eigenvectors.
c
c    Output, double precision D(N), the eigenvalues, in descending order.
c
c    Output, integer IT_NUM, the total number of iterations.
c
c    Output, integer ROT_NUM, the total number of rotations.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision bw(n)
      double precision c
      double precision d(n)
      double precision g
      double precision gapq
      double precision h
      integer i
      integer it_max
      integer it_num
      integer j
      integer k
      integer l
      integer m
      integer p
      integer q
      integer rot_num
      double precision s
      double precision t
      double precision tau
      double precision term
      double precision termp
      double precision termq
      double precision theta
      double precision thresh
      double precision v(n,n)
      double precision w(n)
      double precision zw(n)

      do j = 1, n 
        do i = 1, n 
          v(i,j) = 0.0D+00
        end do
        v(j,j) = 1.0D+00
      end do

      do i = 1, n
        d(i) = a(i,i)
      end do

      do i = 1, n
        bw(i) = d(i)
        zw(i) = 0.0D+00
      end do

      it_num = 0
      rot_num = 0

10    continue

      if ( it_num .lt. it_max ) then

        it_num = it_num + 1
c
c  The convergence threshold is based on the size of the elements in
c  the strict upper triangle of the matrix.
c
        thresh = 0.0D+00
        do j = 1, n
          do i = 1, j - 1
            thresh = thresh + a(i,j) ** 2
          end do
        end do

        thresh = sqrt ( thresh ) / dble ( 4 * n )

        if ( thresh .eq. 0.0D+00 ) then
          go to 20
        end if

        do p = 1, n
          do q = p + 1, n

            gapq = 10.0D+00 * abs ( a(p,q) )
            termp = gapq + abs ( d(p) )
            termq = gapq + abs ( d(q) )
c
c  Annihilate tiny offdiagonal elements.
c
            if ( 4 .lt. it_num .and.
     &           termp .eq. abs ( d(p) ) .and.
     &           termq .eq. abs ( d(q) ) ) then

              a(p,q) = 0.0D+00
c
c  Otherwise, apply a rotation.
c
            else if ( thresh .le. abs ( a(p,q) ) ) then

              h = d(q) - d(p)
              term = abs ( h ) + gapq

              if ( term .eq. abs ( h ) ) then
                t = a(p,q) / h
              else
                theta = 0.5D+00 * h / a(p,q)
                t = 1.0D+00 / 
     &            ( abs ( theta ) + sqrt ( 1.0D+00 + theta * theta ) )
                if ( theta .lt. 0.0D+00 ) then
                  t = - t
                end if
              end if

              c = 1.0D+00 / sqrt ( 1.0D+00 + t * t )
              s = t * c
              tau = s / ( 1.0D+00 + c )
              h = t * a(p,q)
c
c  Accumulate corrections to diagonal elements.
c
              zw(p) = zw(p) - h
              zw(q) = zw(q) + h
              d(p) = d(p) - h
              d(q) = d(q) + h

              a(p,q) = 0.0D+00
c
c  Rotate, using information from the upper triangle of A only.
c
              do j = 1, p - 1
                g = a(j,p)
                h = a(j,q)
                a(j,p) = g - s * ( h + g * tau )
                a(j,q) = h + s * ( g - h * tau )
              end do

              do j = p + 1, q - 1
                g = a(p,j)
                h = a(j,q)
                a(p,j) = g - s * ( h + g * tau )
                a(j,q) = h + s * ( g - h * tau )
              end do

              do j = q + 1, n
                g = a(p,j)
                h = a(q,j)
                a(p,j) = g - s * ( h + g * tau )
                a(q,j) = h + s * ( g - h * tau )
              end do
c
c  Accumulate information in the eigenvector matrix.
c
              do j = 1, n
                g = v(j,p)
                h = v(j,q)
                v(j,p) = g - s * ( h + g * tau )
                v(j,q) = h + s * ( g - h * tau )
              end do

              rot_num = rot_num + 1

            end if

          end do
        end do

        do i = 1, n
          bw(i) = bw(i) + zw(i)
          d(i) = bw(i)
          zw(i) = 0.0D+00
        end do

        go to 10

      end if

20    continue
c
c  Restore upper triangle of input matrix.
c
      do j = 1, n
        do i = 1, j - 1
          a(i,j) = a(j,i)
        end do
      end do
c
c  Ascending sort the eigenvalues and eigenvectors.
c
      do k = 1, n - 1

        m = k

        do l = k + 1, n
          if ( d(l) .lt. d(m) ) then
            m = l
          end if
        end do

        if ( m .ne. k ) then

          t    = d(m)
          d(m) = d(k)
          d(k) = t

          do i = 1, n
            w(i)   = v(i,m)
            v(i,m) = v(i,k)
            v(i,k) = w(i)
          end do

        end if

      end do

      return
      end
