C************************************************************************
        subroutine gtmc(p,npatt,r,patt,mc,nmc,last)
C Finds the column numbers of the missing variables, and stores them
C in the first nmc elements of mc
        integer p,npatt,r(npatt,p),patt,mc(p),nmc,last
        nmc=0
        do 10 j=1,last
           if(r(patt,j).eq.0)then
              nmc=nmc+1
              mc(nmc)=j
           endif
10      continue
        return
        end
C************************************************************************
        subroutine initc(p,c,mc,nmc)
        integer p,c(p),mc(p),nmc
        do 1 j=1,nmc
           c(mc(j))=1
1       continue
        return
        end
C************************************************************************
        subroutine gtdmis(p,d,mc,nmc,dmis)
        integer p,d(p),mc(p),nmc,dmis
        dmis=1
        do 1 j=1,nmc
           dmis=dmis*d(mc(j))
1       continue
        return
        end
C************************************************************************
        subroutine advc(p,c,d,mc,nmc)
C Advances c to next value
        integer p,c(p),d(p),mc(p),nmc
        do 1 j=1,nmc
           if(c(mc(j)).lt.d(mc(j)))then
              c(mc(j))=c(mc(j))+1
              goto 2
           else
              c(mc(j))=1
           endif
1       continue
2       continue
        return
        end
C************************************************************************
        subroutine gtmmis(p,c,mc,nmc,jmp,mmis)
C Calculates mmis from c
        integer p,c(p),mc(p),nmc,jmp(p),mmis
        mmis=0
        do 1 j=1,nmc
           mmis=mmis+(c(mc(j))-1)*jmp(mc(j))
1       continue
        return
        end
C************************************************************************
        subroutine sumc(p,c,mc,nmc,d,jmp,mobs0,dmis,ncells,theta,sum)
        integer p,c(p),mc(p),nmc,d(p),jmp(p),mobs0
        integer m,mmis,dmis,ncells
        double precision sum,theta(ncells)
        call initc(p,c,mc,nmc)
        sum=dfloat(0)
        mmis=0
        do 1 i=1,dmis
           if(i.ne.1)then
              call advc(p,c,d,mc,nmc)
              call gtmmis(p,c,mc,nmc,jmp,mmis)
           endif
           m=mobs0+mmis
           sum=sum+theta(m)
1       continue
        return
        end
C************************************************************************
        subroutine estepc(ncells,theta,x,npatt,p,r,mdpgrp,ngrp,mobs,
     /     nmobs,d,jmp,err,c,mc)
C Performs E-step of EM algorithm. Allocates partially observed data
C to cells of x under specified value of theta.
C If a partially observed unit appears in a cell for which the
C margin of theta is zero, then err is set to one.
C Note: elements of theta need not sum to one.
        integer ncells,npatt,p,r(npatt,p),mdpgrp(npatt),ngrp
        integer mobs(ngrp),nmobs(ngrp),d(p),jmp(p),c(p),mc(p)
        integer patt,grpno,dmis,m,mmis,err,nmc
        double precision sum,theta(ncells),x(ncells)
        do 1 i=1,ncells
           x(i)=dfloat(0)
1       continue
        err=0
        grpno=0
        do 200 patt=1,npatt
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           if(nmc.eq.0)then
              do 50 k=1,mdpgrp(patt)
                 grpno=grpno+1
                 if(theta(mobs(grpno)).le.dfloat(0))then
                    err=1
                    goto 250
                 endif
                 x(mobs(grpno))=x(mobs(grpno))+dfloat(nmobs(grpno))
50            continue
           else
              call gtdmis(p,d,mc,nmc,dmis)
              do 150 k=1,mdpgrp(patt)
                 grpno=grpno+1
                 call sumc(p,c,mc,nmc,d,jmp,mobs(grpno),dmis,
     /                ncells,theta,sum)
                 if(sum.le.dfloat(0))then
                    err=1
                    goto 250
                 endif
                 call initc(p,c,mc,nmc)
                 mmis=0
                 do 130 i=1,dmis
                    if(i.ne.1)then
                       call advc(p,c,d,mc,nmc)
                       call gtmmis(p,c,mc,nmc,jmp,mmis)
                    endif
                    m=mobs(grpno)+mmis
                    x(m)=x(m)+dfloat(nmobs(grpno))*theta(m)/sum
130              continue
150           continue
           endif
200     continue
250     continue
        return
        end
C************************************************************************
        subroutine istepc(ncells,theta,x,npatt,p,r,mdpgrp,ngrp,mobs,
     /     nmobs,d,jmp,err,c,mc)
C Performs I-step of DA algorithm. Allocates partially observed data
C to cells of x at random under specified value of theta.
C If a partially observed unit appears in a cell for which the
C margin of theta is zero, then err is set to one.
C Note: Elements of theta need not sum to one.
        integer ncells,npatt,p,r(npatt,p),mdpgrp(npatt),ngrp
        integer mobs(ngrp),nmobs(ngrp),d(p),jmp(p),c(p),mc(p)
        integer patt,grpno,dmis,m,mmis,err,nmc
        double precision sum,theta(ncells),x(ncells),sum1,u
        do 1 i=1,ncells
           x(i)=dfloat(0)
1       continue
        err=0
        grpno=0
        do 200 patt=1,npatt
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           if(nmc.eq.0)then
              do 50 k=1,mdpgrp(patt)
                 grpno=grpno+1
                 if(theta(mobs(grpno)).le.dfloat(0))then
                    err=1
                    goto 250
                 endif
                 x(mobs(grpno))=x(mobs(grpno))+dfloat(nmobs(grpno))
50            continue
           else
              call gtdmis(p,d,mc,nmc,dmis)
              do 150 k=1,mdpgrp(patt)
                 grpno=grpno+1
                 call sumc(p,c,mc,nmc,d,jmp,mobs(grpno),dmis,
     /                ncells,theta,sum)
                 if(sum.le.dfloat(0))then
                    err=1
                    goto 250
                 endif
C Distribute data at random
                 do 140 l=1,nmobs(grpno)
                    call initc(p,c,mc,nmc)
                    mmis=0
                    u=dble(rangen(0))
                    sum1=0
                    do 130 i=1,dmis
                       if(i.ne.1)then
                          call advc(p,c,d,mc,nmc)
                          call gtmmis(p,c,mc,nmc,jmp,mmis)
                       endif
                       m=mobs(grpno)+mmis
                       sum1=sum1+(theta(m)/sum)
                       if((sum1.gt.u).or.(i.eq.dmis))then
                          x(m)=x(m)+dfloat(1)
                          goto 135
                       endif
130                 continue
135                 continue
140              continue
150           continue
           endif
200     continue
250     continue
        return
        end
C************************************************************************
        subroutine pstep1c(ncells,x,theta,err)
C Performs posterior step of DA. Draws a value theta from
C the Dirichlet posterior distribution of theta, given updated 
C hyperparameters in x. Structural zeroes are 
C flagged by x=-999 and set to zero.
C If any element of x is <= 0 but not -999, then err
C is set to one.
        integer ncells,err
        double precision sum,x(ncells),theta(ncells)
        err=0
        sum=dfloat(0)
        do 5 i=1,ncells
           if(x(i).ne.dfloat(-999))then
              if(x(i).le.dfloat(0))then
                 err=1
                 goto 15
              else
                 theta(i)=dble(gamm(sngl(x(i))))
              endif
              sum=sum+theta(i)
           endif
5       continue
        do 10 i=1,ncells
           if(x(i).ne.dfloat(-999))then
              theta(i)=theta(i)/sum
           else
              theta(i)=dfloat(0)
           endif
10      continue
 15     continue
        return
        end
C************************************************************************
        subroutine llc(ncells,theta,npatt,p,r,mdpgrp,
     /     ngrp,mobs,nmobs,d,jmp,c,mc,ll,err)
C Calculates observed-data loglikelihood function. If a nonzero
C count appears in a cell for which the probability implied by 
C theta is zero, then sets err=1.
        integer ncells,npatt,p,r(npatt,p),mdpgrp(npatt),ngrp
        integer mobs(ngrp),nmobs(ngrp),d(p),jmp(p)
        integer dmis,patt,c(p),mc(p),grpno,nmc,err
        double precision sum,ll,tmp,theta(ncells)
        err=0
        ll=dfloat(0)
        grpno=0
        do 200 patt=1,npatt
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtdmis(p,d,mc,nmc,dmis)
           do 150 k=1,mdpgrp(patt)
              grpno=grpno+1
              call sumc(p,c,mc,nmc,d,jmp,mobs(grpno),dmis,
     /             ncells,theta,sum)
              tmp=dfloat(nmobs(grpno))
              if(sum.le.dfloat(0))then
                 err=1
                 goto 250
              else
                 ll=ll+tmp*log(sum)
              endif
150        continue
200     continue
250     continue
        return
        end
C************************************************************************
        subroutine impc(n,p,x,ncells,theta,npatt,r,mdpgst,mdpgrp,
     /     ngrp,mobs,nmobs,mobsst,d,jmp,c,mc)
C Performs imputation of missing data under theta. Result is
C a completed data matrix x.
        integer n,p,x(n,p),ncells,npatt,mdpgst(npatt)
        integer r(npatt,p),mdpgrp(npatt),ngrp,mobs(ngrp)
        integer nmobs(ngrp),mobsst(ngrp),d(p),jmp(p),mc(p)
        integer dmis,patt,c(p),grpno,mmis,m,drow,dcol,nmc
        double precision sum,sum1,theta(ncells),u
        do 200 patt=1,npatt
           call gtmc(p,npatt,r,patt,mc,nmc,p)
           call gtdmis(p,d,mc,nmc,dmis)
           if(nmc.gt.0)then
C Find the conditional disn of Xmis given Xobs
              do 150 grpno=mdpgst(patt),(mdpgst(patt)+
     /           mdpgrp(patt)-1)
                 call sumc(p,c,mc,nmc,d,jmp,mobs(grpno),dmis,
     /                ncells,theta,sum)
C Distribute data at random
                 drow=mobsst(grpno)-1
                 do 140 l=1,nmobs(grpno)
                    drow=drow+1
                    call initc(p,c,mc,nmc)
                    mmis=0
                    u=dble(rangen(0))
                    sum1=0
                    do 130 i=1,dmis
                       if(i.ne.1)then
                          call advc(p,c,d,mc,nmc)
                          call gtmmis(p,c,mc,nmc,jmp,mmis)
                       endif
                       m=mobs(grpno)+mmis
                       sum1=sum1+(theta(m)/sum)
                       if((sum1.gt.u).or.(i.eq.dmis))then
                          do 95 dcol=1,nmc
                             x(drow,mc(dcol))=c(mc(dcol))
95                        continue  
                          goto 135
                       endif
130                 continue
135                 continue
140              continue
150           continue
           endif   
200     continue
        return
        end
C************************************************************************
        function rangen(init)
        integer a,p,ix,b15,b16,xhi,xalo,leftflo,fhi,k,init
        data a/16807/,b15/32768/,b16/65536/,p/2147483647/
        save ix
        if(init.ne.0) ix=init
        xhi=ix/b16
        xalo=(ix-xhi*b16)*a
        leftflo=xalo/b16
        fhi=xhi*a+leftflo
        k=fhi/b15
        ix=(((xalo-leftflo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0)ix=ix+p
        rangen=float(ix)*4.656612875E-10
        return
        end
C***********************************************************************
        subroutine rngs(seed)
C initializes rangen with seed
        integer seed
        tmp=rangen(seed)
        return
        end
C***********************************************************************
        function gamm(a)
C Generates a random gamma(a) variate. If a>=1, uses the method of 
C Fishman (1976); if 0<a<1, the method of Ahrens (1974)
        real a,u,y,q,e,b,p,u1,lq
        data e/2.718282/
        if(a.ge.1)then
1          continue
           u=rangen(0)
           y=-log(rangen(0))
C           q=(y/exp(y-1))**(a-1)
           lq=(a-1.)*(log(y)-(y-1.))
           q=exp(lq)
           if(u.le.q)then
              gamm=a*y
           else
              goto 1
           endif
        else
2          continue
           u=rangen(0)
           b=(e+a)/e
           p=b*u
           if(p.gt.1) goto 4
3          continue
           x=p**(1/a)
           u1=rangen(0)
           if(u1.gt.(e**(-x)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
4          continue
           x=-log((b-p)/a)
           u1=rangen(0)
           if(u1.gt.(x**(a-1)))then
              goto 2
           else
              gamm=x
              goto 10
           endif
        endif
10      continue
        return
        end
C***********************************************************************
        subroutine mgstepc(ncells,oldtheta,newtheta,npatt,p,r,mdpgrp,
     /     ngrp,mobs,nmobs,d,jmp,c,mc,dtable,prior,err,sj)
C Performs one step of monotone DA. Randomly allocates partially 
C observed data to cells of newtheta, assuming that theta=oldtheta.
C A structural zero is indicated by oldtheta=0 and prior=-999.
C If any observation is allocated to a structural or random zero in the
C I-step, sets err to one. This indicates that the observed data are
C inconsistent with the starting value and/or prior.
C If any prior count + cell count < 0, sets err to one.
C If posterior is improper, sets err to 2.
        integer ncells,npatt,p,r(npatt,p),mdpgrp(npatt),ngrp
        integer mobs(ngrp),nmobs(ngrp),d(p),jmp(p)
        integer dmis,patt,c(p),mc(p),grpno,mmis,m,nmc,lp,fp
        integer outer,inner,ii,jmpii,nii,mii,err,sj(p),zflag
        double precision oldtheta(ncells),newtheta(ncells)
        double precision prior(ncells),u,sum,sum1,sum2,dtable(ncells)
        err=0
C put prior counts in dtable, initialize newtheta to all ones
        do 5 i=1,ncells
           dtable(i)=prior(i)
           newtheta(i)=dfloat(1)
5       continue
        grpno=0
        do 400 j=p,1,-1
           if(j.eq.p) then
              fp=1
           else
              fp=sj(j+1)+1
           endif
           lp=sj(j)
C do the I-step
           do 200 patt=fp,lp
              call gtmc(p,npatt,r,patt,mc,nmc,j)
              call gtdmis(p,d,mc,nmc,dmis)
C Find the conditional disn of Ymis* given Yobs
              do 150 k=1,mdpgrp(patt)
                 grpno=grpno+1
                 call sumc(p,c,mc,nmc,d,jmp,mobs(grpno),dmis,
     /             ncells,oldtheta,sum)
C Distribute data at random by table sampling
                 do 140 l=1,nmobs(grpno)
                    call initc(p,c,mc,nmc)
                    mmis=0
                    u=dble(rangen(0))
                    sum1=dfloat(0)
                    do 130 i=1,dmis
                       if(i.ne.1)then
                          call advc(p,c,d,mc,nmc)
                          call gtmmis(p,c,mc,nmc,jmp,mmis)
                       endif
                       m=mobs(grpno)+mmis
                       sum1=sum1+(oldtheta(m)/sum)
                       if((sum1.gt.u).or.(i.eq.dmis))then
                          if(oldtheta(m).eq.dfloat(0))then
                             err=1
                             goto 450
                          endif
                          dtable(m)=dtable(m)+dfloat(1)
                          goto 135
                       endif
130                 continue
135                 continue
140              continue
150           continue
200        continue
C Now draw from the jth posterior factor, while collapsing dtable
C and oldtheta
           jmpii=jmp(j)*d(j)
           nii=ncells/jmpii
           do 350 outer=1,jmp(j)
              m=outer-jmp(j)
              sum=dfloat(0)
              sum1=dfloat(0)
              sum2=dfloat(0)
              zflag=1
              do 250 inner=1,d(j)
                 m=m+jmp(j)
                 if(dtable(m).eq.dfloat(-999))then
                    u=dfloat(0)
                 else
                    if(dtable(m).le.dfloat(0))then
                       err=2
                       goto 450
                    else
                       u=dble(gamm(sngl(dtable(m))))
                       zflag=0
                       sum1=sum1+dtable(m)
                    endif
                 endif
                 sum=sum+u
                 sum2=sum2+oldtheta(m)
                 mii=m-jmpii
                 do 230 ii=1,nii
                    mii=mii+jmpii
                    newtheta(mii)=newtheta(mii)*u
230              continue
250           continue
              oldtheta(outer)=sum2
              if(zflag.eq.1)then
                 dtable(outer)=dfloat(-999)
              else
                 dtable(outer)=sum1
              endif
              m=outer-jmp(j)
              do 300 inner=1,d(j)
                 m=m+jmp(j)
                 dtable(m)=sum1
                 mii=m-jmpii
                 do 280 ii=1,nii
                    mii=mii+jmpii
                    newtheta(mii)=newtheta(mii)/sum
280              continue
300           continue
350        continue
400     continue
450     continue
        return
        end
C************************************************************************
        subroutine gtntab(ncon,con,ntab)
C find number of marginal tables to be fit
        integer ncon,con(ncon),ntab,flag,posn
        ntab=0
        flag=0
        do 10 posn=1,ncon
           if((con(posn).ne.0).and.(flag.eq.0))flag=1
           if((con(posn).eq.0).and.(flag.eq.1))then
              ntab=ntab+1
              flag=0
           endif
           if((flag.eq.1).and.(posn.eq.ncon))ntab=ntab+1
10      continue
        return
        end
C************************************************************************
        subroutine gtmarg(ncon,con,posn,p,marg,nmarg)
C extract the next set of margins to be fit, store them in the first 
C nmarg elements of marg. For first set, posn should be 0.
        integer ncon,con(ncon),p,marg(p),nmarg,posn
1       continue
           posn=posn+1
           if(con(posn).eq.0)goto 1
        nmarg=0
11      continue
          if(con(posn).eq.0)goto 12
          nmarg=nmarg+1
          marg(nmarg)=con(posn)
          if(posn.eq.ncon)goto 12
          posn=posn+1
          goto 11
12      continue
        return
        end
C************************************************************************
        subroutine gtrest(p,marg,nmarg,rest,nrest)
        integer p,marg(p),nmarg,rest(p),nrest,flag
        nrest=0
        do 20 j=1,p
           flag=0
           do 10 k=1,nmarg
              if(j.eq.marg(k))then
                 flag=1
                 goto 15
              endif
10         continue
15         continue
           if(flag.eq.0)then
              nrest=nrest+1
              rest(nrest)=j
           endif
20      continue
        return
        end
C************************************************************************
        subroutine ipf(ncells,table,fit,ncon,con,p,d,jmp,c,marg,rest,
     /     eps)
C Performs one step of iterative proportional fitting, raking fit
C to the margins of table. Margins to be fit are specified in con.
        integer ncells,ncon,con(ncon),p,d(p),jmp(p),c(p),marg(p)
        integer rest(p),nmarg,nrest,m,mmarg,mrest,posn,ntab,tabno
        integer out,in,dmarg,drest
        double precision sumt,sumf,table(ncells),fit(ncells),eps
        call gtntab(ncon,con,ntab)
        posn=0
        do 100 tabno=1,ntab
           call gtmarg(ncon,con,posn,p,marg,nmarg)
           call gtrest(p,marg,nmarg,rest,nrest)
           call gtdmis(p,d,marg,nmarg,dmarg)
           drest=ncells/dmarg
           call initc(p,c,marg,nmarg)
           mmarg=1
           do 90 out=1,dmarg
              if(out.ne.1)then
                 call advc(p,c,d,marg,nmarg)
                 call gtmmis(p,c,marg,nmarg,jmp,mmarg)
                 mmarg=mmarg+1
              endif
              call sum2c(p,c,rest,nrest,d,jmp,mmarg,drest,ncells,
     /             table,sumt,fit,sumf)
              call initc(p,c,rest,nrest)
              if(sumf.ne.0)then
                 mrest=0
                 do 80 in=1,drest   
                    if(in.ne.1)then
                       call advc(p,c,d,rest,nrest)
                       call gtmmis(p,c,rest,nrest,jmp,mrest)
                    endif
                    m=mmarg+mrest
                    if(fit(m).ge.eps) then
                       fit(m)=fit(m)*(sumt/sumf)
                    else
                       fit(m)=dfloat(0)
                    endif
80               continue
              endif
90         continue
100     continue
        return
        end
C************************************************************************
        subroutine bipf(ncells,table,theta,prior,ncon,con,p,d,jmp,c,
     /     marg,rest,err)
C Performs one cycle of Bayesian ipf. Cell counts are in table
C and starting value in theta.  Prior hyperparameters are in prior,
C with structural zeros denoted by -999. Replaces theta with an
C updated value. 
        integer ncells,ncon,con(ncon),p,d(p),jmp(p),c(p),marg(p)
        integer rest(p),nmarg,nrest,m,mmarg,mrest,posn,ntab,tabno
        integer out,in,dmarg,drest,err,zflag
        double precision sumt,sumf,table(ncells),theta(ncells),g
        double precision prior(ncells),sum3
        call gtntab(ncon,con,ntab)
        err=0
        posn=0
        do 100 tabno=1,ntab
           sum3=dfloat(0)
           call gtmarg(ncon,con,posn,p,marg,nmarg)
           call gtrest(p,marg,nmarg,rest,nrest)
           call gtdmis(p,d,marg,nmarg,dmarg)
           drest=ncells/dmarg
           call initc(p,c,marg,nmarg)
           mmarg=1
           do 90 out=1,dmarg
              if(out.ne.1)then
                 call advc(p,c,d,marg,nmarg)
                 call gtmmis(p,c,marg,nmarg,jmp,mmarg)
                 mmarg=mmarg+1
              endif
              zflag=0
              call sum3c(p,c,rest,nrest,d,jmp,mmarg,drest,ncells,
     /             table,sumt,theta,sumf,prior,zflag)
              call initc(p,c,rest,nrest)
              if(sumt.le.0)then
                 err=1
                 goto 150
              endif
              if(zflag.eq.1)then
                 g=dble(gamm(sngl(sumt)))+dble(1e-20)
                 sum3=sum3+g
              endif
              mrest=0
              do 80 in=1,drest   
                 if(in.ne.1)then
                    call advc(p,c,d,rest,nrest)
                    call gtmmis(p,c,rest,nrest,jmp,mrest)
                 endif
                 m=mmarg+mrest
                 theta(m)=theta(m)*g/sumf
 80           continue
 90        continue
           do 95 m=1,ncells
              theta(m)=theta(m)/sum3
 95        continue
100     continue
150     continue
        return
        end
C************************************************************************
        subroutine sum2c(p,c,mc,nmc,d,jmp,mobs,dmis,
     /      ncells,table1,sum1,table2,sum2)
        integer p,c(p),mc(p),nmc,d(p),jmp(p),mobs
        integer m,mmis,dmis,ncells
        double precision sum1,sum2,table1(ncells),table2(ncells)
        call initc(p,c,mc,nmc)
        sum1=dfloat(0)
        sum2=dfloat(0)
        mmis=0
        do 1 i=1,dmis
           if(i.ne.1)then
              call advc(p,c,d,mc,nmc)
              call gtmmis(p,c,mc,nmc,jmp,mmis)
           endif
           m=mobs+mmis
           sum1=sum1+table1(m)
           sum2=sum2+table2(m)
1       continue
        return
        end
C************************************************************************
        subroutine sum3c(p,c,mc,nmc,d,jmp,mobs,dmis,
     /      ncells,table1,sum1,table2,sum2,prior,zflag)
        integer p,c(p),mc(p),nmc,d(p),jmp(p),mobs
        integer m,mmis,dmis,ncells,zflag
        double precision sum1,sum2,table1(ncells),table2(ncells),
     /    prior(ncells)
        call initc(p,c,mc,nmc)
        sum1=dfloat(0)
        sum2=dfloat(0)
        mmis=0
        do 1 i=1,dmis
           if(i.ne.1)then
              call advc(p,c,d,mc,nmc)
              call gtmmis(p,c,mc,nmc,jmp,mmis)
           endif
           m=mobs+mmis
           sum2=sum2+table2(m)
           if(prior(m).ne.dfloat(-999)) then
              sum1=sum1+table1(m)+prior(m)
              zflag=1
           endif
1       continue
        return
        end
C************************************************************************
        subroutine sjn(p,npatt,r,sj)
C computes s_j, the number of the last missingness pattern for which
C the jth variable needs to be imputed to complete the monotone pattern
        integer p,npatt,r(npatt,p),sj(p),patt,tmp
        do 30 j=1,p
           patt=npatt+1
10         continue
           patt=patt-1
           if(patt.ge.1)then
              if(r(patt,j).eq.0)goto 10
           endif
           sj(j)=patt
30      continue
        tmp=sj(p)
        do 40 j=p-1,1,-1
           sj(j)=max0(sj(j),tmp)
           tmp=sj(j)
40      continue
        return
        end
C***********************************************************************
