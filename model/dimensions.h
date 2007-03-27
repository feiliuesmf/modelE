c-----------------------------------------------------------------------------
      integer idm,jdm,kdm,ms,iia,jja,iio,jjo,ntrcr,iold
      real equato
      parameter (idm=195,jdm=180,kdm=20,ms=15,ntrcr=1,equato=115.
     .          ,iold=181)
      parameter (iia=72,jja=46,iio=idm,jjo=jdm)
c
c --- ms-1  = max. number of interruptions of any grid row or column by land
c
      integer ii,jj,kk,ii1,nlato,nlongo,nlatn,nlongn
      parameter (ii=idm,jj=jdm,kk=kdm,ii1=ii-1)
      parameter (nlato=1,nlongo=1,nlatn=idm,nlongn=jdm)
c
c --- information in common block gindex keeps do loops from running into land
      common/gindex/ip(idm,jdm),iu(idm,jdm),iv(idm,jdm),iq(idm,jdm),
     .ifp(jdm,ms),ilp(jdm,ms),isp(jdm),jfp(idm,ms),jlp(idm,ms),jsp(idm),
     .ifq(jdm,ms),ilq(jdm,ms),isq(jdm),jfq(idm,ms),jlq(idm,ms),jsq(idm),
     .ifu(jdm,ms),ilu(jdm,ms),isu(jdm),jfu(idm,ms),jlu(idm,ms),jsu(idm),
     .ifv(jdm,ms),ilv(jdm,ms),isv(jdm),jfv(idm,ms),jlv(idm,ms),jsv(idm),
     .msk(idm,jdm)
c
      integer ip,iu,iv,iq,
     .        ifp,ilp,isp,jfp,jlp,jsp,ifq,ilq,isq,jfq,jlq,jsq,
     .        ifu,ilu,isu,jfu,jlu,jsu,ifv,ilv,isv,jfv,jlv,jsv,msk
      integer lp
c
      common /linepr/ lp
c-----------------------------------------------------------------------------
