      /*
    if(recon->iconj[i] == 0)
    {
      numerx = conj(recon->gx_ft[recon->ifs[i]]) * recon->sx_ft[i];
      numery = conj(recon->gy_ft[recon->ifs[i]]) * recon->sy_ft[i];
      numer = numerx + numery;
      est = numer * recon->gd_ft[recon->ifs[i]];
      recon->est_ft[i] = est;

    }else{
      numerx = conj(recon->gx_ft[recon->ifs[i]]) * recon->sx_ft[i];
      numery = conj(recon->gy_ft[recon->ifs[i]]) * recon->sy_ft[i];
      numer = numerx + numery;
      est = numer * recon->gd_ft[recon->ifs[i]];
      recon->est_ft[i] = est;
      // recon->est_ft[i] = (conj(recon->gx_ft[recon->ifs[i]]) * conj(recon->sx_ft[i])
//         + conj(recon->gy_ft[recon->ifs[i]]) * conj(recon->sy_ft[i]))
//         * recon->gd_ft[recon->ifs[i]];
    }
    */
      /*
    if(cabs(recon->est_ft[i]) > 0.001)
    {
      x = i % recon->nf;
      y = i / recon->nf;
      printf("c:  ");
      printf("sx[%d,%d]=%g%+gi ",y,x,creal(recon->sx_ft[i]), cimag(recon->sx_ft[i]));
      printf("sy[%d,%d]=%g%+gi ",y,x,creal(recon->sy_ft[i]), cimag(recon->sy_ft[i]));
      printf("es[%d,%d]=%g%+gi ",y,x,creal(recon->est_ft[i]), cimag(recon->est_ft[i]));
      printf("\n");
      printf("c:  ");
      printf("nx[%d,%d]=%g%+gi ",y,x, creal(numerx), cimag(numerx));
      printf("ny[%d,%d]=%g%+gi ",y,x, creal(numery), cimag(numery));
      printf("ny[%d,%d]=%g%+gi ",y,x, creal(numer), cimag(numer));
      printf("est[%d,%d]=%g%+gi ",y,x, creal(est), cimag(est));
      printf("\n");
      printf("c:  ");
      printf("gx[%d,%d]=%g%+gi ",y,x,creal(recon->gx_ft[recon->ifs[i]]), cimag(recon->gx_ft[recon->ifs[i]]));
      printf("gy[%d,%d]=%g%+gi ",y,x,creal(recon->gy_ft[recon->ifs[i]]), cimag(recon->gy_ft[recon->ifs[i]]));
      printf("gd[%d,%d]=%g%+gi ",y,x,creal(recon->gd_ft[recon->ifs[i]]), cimag(recon->gd_ft[recon->ifs[i]]));
      printf("\n");
    }
      */