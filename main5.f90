!Program: BEM_mesh_quad_dis
!Dveloped by Caio Moura
!UNICAMP, 2021
!Function: It is to build a discontinuous quadratic boundary element mesh,
!apply the boundary conditions on these elements and solve the BEM Potencial problem

PROGRAM BEM_mesh_quad_dis

    IMPLICIT NONE
    INTEGER :: count_0, count_1, count_rate, count_max, i, j, k, n, m, sumk, num_arestas_totais, &
    num_arestas_cc_temp, num_arestas_cc_fluxo, ne, s, npg, npglog, xstepintern, ystepintern, &
    neinternal
    DOUBLEPRECISION :: xo, xf, yo, yf, xee, yee, incrementoy, incrementox, xd, yd, &
    kmaterial, dFF1, dFF2, dFF3, x1, x2, x3, y1, y2, y3, xpi, ypi, dxdcsi, &
    dydcsi, jac, nx, ny, Rx, Ry, R, Tsolfun, Qsolfun, csino, &
    newcsilog, pi, RB, RB1, RB2, time_init, time_final, elapsed_time, &
    btemp, somanorma, normau, prodinter, somats, newcsi, Jcsi, &
    fx, fy, ff1, ff2, ff3, k1, k2, k3, intG1, intG2, intG3, intH1, &
    intH2, intH3, intnsg1A, intnsg2A, intnsg3A, intsg1A, intsg2A, intsg3A, &
    intnsg1B, intnsg2B, intnsg3B, intsg1B, intsg2B, intsg3B, xinc, yinc, &
    xintersum, yintersum, dx, dy
    DOUBLEPRECISION, ALLOCATABLE :: num_vertice(:), x_vertice(:), y_vertice(:), &
    valor_temp_aresta(:), valor_fluxo_aresta(:), xe(:), ye(:), xe_temp(:), ye_temp(:), &
    xno(:), yno(:), xnofinal(:), ynofinal(:), csi(:), weights(:), termoconhe(:), &
    csilog(:), weightslog(:), FF(:), G_temp(:,:), H_temp(:,:), G(:,:), H(:,:),&
    termoconhearesta(:), inversorvector1(:), inversorvector2(:), x(:), b(:), &
    A(:,:), coluna_atual(:), Q(:,:), QT(:,:), proatual(:), u(:), y(:), Rsystem(:,:), &
    xelefinal(:), yelefinal(:), NF(:), inode(:), mnode(:), fnode(:), Qbound(:), Tbound(:), &
    Hd(:,:), Hdx(:,:), Hdy(:,:), Gd(:,:), Gdx(:,:), Gdy(:,:), xinternal(:), yinternal(:), &
    Tin(:),  Qinx(:),  Qiny(:)
    INTEGER, ALLOCATABLE :: xglobal(:), yglobal(:), ref_aresta(:), aresta_cc_temp(:),&
    num_aresta(:), num_vertice_inicial(:), num_vertice_final(:), aresta_cc_fluxo(:), &
    num_ele_aresta(:)

    !Starting timer
    CALL system_clock(count_0, count_rate, count_max)
    time_init = count_0*1.0/count_rate

    PRINT*," "
    PRINT*,"Mesh_quad_dis was started."
    PRINT*," "
    PRINT*,"TXT file reader."
    !Reader of key information 
    PRINT*," -Reading .txt file."
    OPEN(unit=100,file='entrada.txt',status='old')
    READ(100,FMT="(/////I20.7)") num_arestas_totais
    READ(100,FMT='(//I20.7)') num_arestas_cc_temp
    READ(100,FMT='(//I20.7)') num_arestas_cc_fluxo
    READ(100,FMT='(//)') 

    !Allocate in function of "n"
    n=num_arestas_totais
    ALLOCATE (num_vertice(1:n), x_vertice(1:n), y_vertice(1:n), num_aresta(1:n), &
    num_vertice_inicial(1:n), num_vertice_final(1:n), ref_aresta(1:n), num_ele_aresta(1:n),&
    aresta_cc_temp(1:n), valor_temp_aresta(1:n), aresta_cc_fluxo(1:n), valor_fluxo_aresta(1:n),&
    xglobal(1:n), yglobal(1:n), csi(1:4), weights(1:4), csilog(1:10), weightslog(1:10), &
    FF(1:3), inversorvector1(1:n), NF(1:3))
                                               
    !Boundary conditions reader 
    DO i=1,num_arestas_totais !Vertices reader 
        READ(100,FMT='(F2.0,1x,F2.0,1x,F2.0)')num_vertice(i), x_vertice(i), y_vertice(i)
    END DO
    READ(100,FMT='(//)')
    DO i=1,num_arestas_totais !Vertices connectivity reader
        READ(100,FMT='(I2,1x,I2,1x,I2)')num_aresta(i), num_vertice_inicial(i), num_vertice_final(i)
    END DO
    READ(100,FMT='(//)')
    DO i=1,num_arestas_totais !Mesh refinement reader 
        READ(100,FMT='(I2,1x,I2)')ref_aresta(i), num_ele_aresta(i)
    END DO
    READ(100,FMT='(//)')
    DO i=1,num_arestas_cc_temp !BC temperature reader 
        READ(100,FMT='(I2,1x,F6.0)')aresta_cc_temp(i), valor_temp_aresta(i)
    END DO
    READ(100,FMT='(//)')
    DO i=1,num_arestas_cc_fluxo !BC heat flux reader 
        READ(100,FMT='(I2,1x,F3.0)')aresta_cc_fluxo(i), valor_fluxo_aresta(i)
    END DO
    READ(100,FMT='(//F10.0)')kmaterial
    READ(100,FMT='(///I2,1x,I2)')xstepintern, ystepintern
    neinternal=xstepintern*ystepintern
    
    !Number of elements calculator 
    ne=0;
    DO i=1,n 
        ne=ne+num_ele_aresta(i)
    END DO
    ALLOCATE (xe(1:ne),ye(1:ne), xe_temp(1:ne), ye_temp(1:ne), xno(1:ne*3), yno(1:ne*3), &
    xnofinal(1:ne*3), ynofinal(1:ne*3), G(1:ne*3,1:ne*3), H(1:ne*3,1:ne*3), &
    G_temp(1:ne*3,1:ne*3), H_temp(1:ne*3,1:ne*3), termoconhearesta(1:ne*3), &
    termoconhe(1:ne*3), inversorvector2(1:ne*3), x(1:ne*3), b(1:ne*3), &
    coluna_atual(1:ne*3), A(1:ne*3,1:ne*3), Q(1:ne*3,1:ne*3), proatual(1:ne*3), &
    u(1:ne*3), QT(1:ne*3,1:ne*3), y(1:ne*3), Rsystem(1:ne*3,1:ne*3), &
    xelefinal(1:ne*3), yelefinal(1:ne*3), inode(1:ne), mnode(1:ne), fnode(1:ne), &
    Qbound(1:ne*3), Tbound(1:ne*3), xinternal(1:neinternal), yinternal(1:neinternal),&
    Hd(1:neinternal,1:ne*3), Hdx(1:neinternal,1:ne*3), Hdy(1:neinternal,1:ne*3), &
    Gd(1:neinternal,1:ne*3), Gdx(1:neinternal,1:ne*3), Gdy(1:neinternal,1:ne*3), &
    Tin(1:neinternal), Qinx(1:neinternal), Qiny(1:neinternal))
    PRINT*, " -Number of elements:",ne,"elements"
    
    !Mesh building
    k=1; yee=0; xee=0;
    DO i=1,n !Building the elements on edges
        !Defining X coordinates
        xo=x_vertice(num_vertice_inicial(i))
        xf=x_vertice(num_vertice_final(i))
        !Defining Y coordinates
        yo=y_vertice(num_vertice_inicial(i))
        yf=y_vertice(num_vertice_final(i))
        incrementox=(xf-xo)/num_ele_aresta(i) 
        incrementoy=(yf-yo)/num_ele_aresta(i)
        !The lasts values stored on vectors xe(:) e y(:) are the origin coordinates
        m=k-1+num_ele_aresta(i)
        sumk=0;
        DO j=k,m
            xee=xee+incrementox
            yee=yee+incrementoy
            xe_temp(j)=xee
            ye_temp(j)=yee
            sumk=sumk+1
        END DO
        k=k+sumk
        !Positioning the origin on the start of xe(:) and ye(:) computational vectors 
        DO j=1,ne 
            IF (j.EQ.1) THEN
                xe(j)=xe_temp(ne)
                ye(j)=ye_temp(ne)
            ELSE
                xe(j)=xe_temp(j-1)
                ye(j)=ye_temp(j-1)
            END IF
        END DO
    END DO 
    m=1;

    !Calculating the X and Y nodes coordinates 
    DO i=1,ne*3,3
        IF (i.LT.4) THEN
            !A
            xno(i)=xe(ne)+(((xe(1)-xe(ne))/2)*(1.0/3.0))
            yno(i)=ye(ne)+(((ye(1)-ye(ne))/2)*(1.0/3.0))
            !B
            xno(i+1)=xno(i)+((xe(1)-xe(ne))/2)*(2.0/3.0)
            yno(i+1)=yno(i)+((ye(1)-ye(ne))/2)*(2.0/3.0)
            !C
            xno(i+2)=xno(i+1)+((xe(1)-xe(ne))/2)*(2.0/3.0)
            yno(i+2)=yno(i+1)+((ye(1)-ye(ne))/2)*(2.0/3.0)
        ELSEIF (i.GT.1) THEN
            m=((i-1)/3)+1
            k=m-1
            !A
            xno(i)=xe(k)+(((xe(m)-xe(k))/2)*(1.0/3.0))
            yno(i)=ye(k)+(((ye(m)-ye(k))/2)*(1.0/3.0))
            !B
            xno(i+1)=xno(i)+((xe(m)-xe(k))/2)*(2.0/3.0)
            yno(i+1)=yno(i)+((ye(m)-ye(k))/2)*(2.0/3.0)
            !C
            xno(i+2)=xno(i+1)+((xe(m)-xe(k))/2)*(2.0/3.0)
            yno(i+2)=yno(i+1)+((ye(m)-ye(k))/2)*(2.0/3.0)
        END IF
    END DO
    DO i=4,ne*3
        xnofinal(i-3)=xno(i)
        ynofinal(i-3)=yno(i)
    END DO 
    DO i=1,3
        xnofinal((ne*3)-3+i)=xno(i)
        ynofinal((ne*3)-3+i)=yno(i)
    END DO 
    s=1;
    DO i=1,ne !To define the node positions 
        inode(i)=s
        mnode(i)=s+1
        fnode(i)=s+2
        s=s+3;
    END DO

    !Defining the elements points X and Y coordinates to use in rA and rB
    k=1;
    DO i=1,ne  
        !A
        xelefinal(k)=xe(i)
        yelefinal(k)=ye(i)
        !B
        xelefinal(k+1)=xelefinal(k)+((xe(i+1)-xe(i))/2)
        yelefinal(k+1)=yelefinal(k)+((ye(i+1)-ye(i))/2)
        !C
        xelefinal(k+2)=xelefinal(k+1)+((xe(i+1)-xe(i))/2)
        yelefinal(k+2)=yelefinal(k+1)+((ye(i+1)-ye(i))/2)
        k=k+3
    END DO

    !Building internal points to a square or rectangle 
    s=1;
    xinc=(x_vertice(s+1)-x_vertice(s))/(xstepintern+1)
    yinc=(y_vertice(s+2)-y_vertice(s+1))/(ystepintern+1)
    xinternal(:)=0; yinternal(:)=0; xintersum=0; yintersum=0;
    s=0;
    DO i=1,ystepintern
        yintersum=yintersum+yinc
        xintersum=0;
        DO j=1,xstepintern
            s=s+1
            xintersum=xintersum+xinc
            xinternal(s)=xintersum
            yinternal(s)=yintersum
        END DO
    END DO
 
    print*,"INC X", xinc
    PRINT*,"INC Y", yinc
    print*,"Xin",xinternal(:)
    print*,""
    print*,"Yin",yinternal(:)

    !BEA formulation
    PRINT*," "
    PRINT*,"BEM Formulation to mesh_quad_dis."
    PRINT*, " -Building [H] and [G]."
    !Gauss points and weights (regular)
    !csi(1)=-0.3399810435848563; csi(2)=0.3399810435848563; csi(3)=-0.8611363115940526;
    !csi(4)=0.8611363115940526; weights(1)=0.6521451548625461; weights(2)=0.6521451548625461; 
    !weights(3)=0.3478548451374538; weights(4)=0.3478548451374538;
    !npg=4;
    csi(1)=-0.932469514203152;
    csi(2)=0.932469514203152;
    csi(3)=-0.661209386466265;
    csi(4)=0.661209386466265;
    csi(5)=-0.238619186083197;
    csi(6)=0.238619186083197;
    weights(1)=0.171324492379107;
    weights(2)=0.171324492379107;
    weights(3)=0.360761573048139;
    weights(4)=0.360761573048139;
    weights(5)=0.467913934572691
    weights(6)=0.467913934572691
    npg=6

    !Gauss points and weights (logarithmic singularity)
    csilog(1)=0.0414484801993832; csilog(2)=0.245474914320602; csilog(3)=0.556165453560278;
    csilog(4)=0.848982394532985; weightslog(1)=0.383464068145135; weightslog(2)=0.386875317774762; 
    weightslog(3)=0.190435126950142; weightslog(4)=0.0392254871299598; 
    npglog=4;
    ![G]=[H]=[0] to start the solution
    G_temp(:,:)=0;
    H_temp(:,:)=0;
    pi=3.14159265358979;
    DO i=1,ne*3 !To define the source point 
        xd=xnofinal(i)
        yd=ynofinal(i)
        s=1;
        DO j=1,ne !To define the element
            x1=xelefinal(s); x2=xelefinal(s+1); x3=xelefinal(s+2)
            y1=yelefinal(s); y2=yelefinal(s+1); y3=yelefinal(s+2)
            IF ((i.NE.inode(j)).AND.(i.NE.mnode(j)).AND.(i.NE.fnode(j))) THEN !Regular integration
                intG1=0; intG2=0; intG3=0; intH1=0; intH2=0; intH3=0; 
                DO m=1,npg !Gauss regular integration
                    FF1=(csi(m)**2-csi(m))/2.      
                    FF2=1.-csi(m)**2               
                    FF3=(csi(m)**2+csi(m))/2.      
                    k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2); !Node shape function 
                    NF(1)=(csi(m)-k2)*csi(m)*k2/k3
                    NF(2)=(k1+k2)*(k1*k2+(k2-k1)*csi(m)-csi(m)**2)/k3      
                    NF(3)=(k1+csi(m))*csi(m)*k1/k3   
                    dFF1=(2.*csi(m)-1)/2.  
                    dFF2=-2.*csi(m)         
                    dFF3=(2.*csi(m)+1)/2.  
                    xpi=(FF1*x1)+(FF2*x2)+(FF3*x3);            
                    ypi=(FF1*y1)+(FF2*y2)+(FF3*y3);             
                    Rx=(xpi-xd)
                    Ry=(ypi-yd)
                    R=((Rx**2)+(Ry**2))**0.5              !Distance between source point and field point (Radius) 
                    fx=((x1-2.*x2+x3)*csi(m)+(x3-x1)/2.); 
                    fy=((y1-2.*y2+y3)*csi(m)+(y3-y1)/2.);
                    jac=((fx**2)+(fy**2))**0.5;           !Jacobian
                    nx=fy/jac;                            !X component of the normal vector 
                    ny=-fx/jac;                           !Y component of the normal vector                                 
                    Tsolfun=(-1./(2.*pi*kmaterial))*log(R);                    !Temperature solfun
                    Qsolfun=(1./(2.*pi*kmaterial))*((Rx*nx)+(Ry*ny))/(R**2);   !Flux solfun
                    intG1=intG1+((Tsolfun*(NF(1))*jac*weights(m)));
                    intG2=intG2+((Tsolfun*(NF(2))*jac*weights(m)));
                    intG3=intG3+((Tsolfun*(NF(3))*jac*weights(m)));
                    intH1=intH1+((Qsolfun*(NF(1))*jac*weights(m)));
                    intH2=intH2+((Qsolfun*(NF(2))*jac*weights(m)));
                    intH3=intH3+((Qsolfun*(NF(3))*jac*weights(m)));
                END DO
                H_temp(i,s)=intH1
                H_temp(i,s+1)=intH2
                H_temp(i,s+2)=intH3
                G_temp(i,s)=intG1
                G_temp(i,s+1)=intG2
                G_temp(i,s+2)=intG3
            ELSE !Singularity integration
                IF (i.EQ.inode(j)) THEN
                    csino=-2.0/3.0
                ELSEIF (i.EQ.mnode(j)) THEN
                    csino=0
                ELSE IF (i.EQ.fnode(j)) THEN
                    csino=2.0/3.0    
                END IF
                intnsg1A=0; intnsg2A=0; intnsg3A=0; !Regular Gauss integration A - After xd
                DO m=1,npg !-1 to 1 
                    Jcsi=(1.-csino)/2.                                  !Qsi jacobian 
                    newcsi=Jcsi*csi(m)+(1.+csino)/2.                    !Qsi(eta)
                    RB1=((x3-x1)+(newcsi+csino)*(x1-2.*x2+x3))**2
                    RB2=((y3-y1)+(newcsi+csino)*(y1-2.*y2+y3))**2
                    RB=(RB1+RB2)**0.5;                                  !RB 
                    k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2); 
                    NF(1)=(newcsi-k2)*newcsi*k2/k3                      !Shape function of functional node 1 
                    NF(2)=(k1+k2)*(k1*k2+(k2-k1)*newcsi-newcsi**2)/k3   !Shape function of functional node 2
                    NF(3)=(k1+newcsi)*newcsi*k1/k3                      !Shape function of functional node 3
                    dxdcsi=((x1-2.*x2+x3)*newcsi+(x3-x1)/2.)
                    dydcsi=((y1-2.*y2+y3)*newcsi+(y3-y1)/2.)
                    jac=((dxdcsi**2)+(dydcsi**2))**0.5;                 !Jacobian
                    intnsg1A=intnsg1A+(log(RB)*NF(1)*jac*Jcsi*weights(m))
                    intnsg2A=intnsg2A+(log(RB)*NF(2)*jac*Jcsi*weights(m))
                    intnsg3A=intnsg3A+(log(RB)*NF(3)*jac*Jcsi*weights(m))
                END DO
                intsg1A=0; intsg2A=0; intsg3A=0; !Logarithmic Gauss integration A - After xd
                DO m=1,npglog !0 to 1
                    Jcsi=(1.-csino);                                        !Qsi jacobian
                    newcsilog=Jcsi*csilog(m)+csino                          !Qsi(eta)
                    k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2);  
                    NF(1)=(newcsilog-k2)*newcsilog*k2/k3                    !Shape function of functional node 1 
                    NF(2)=(k1+k2)*(k1*k2+(k2-k1)*newcsilog-newcsilog**2)/k3 !Shape function of functional node 2    
                    NF(3)=(k1+newcsilog)*newcsilog*k1/k3                    !Shape function of functional node 3
                    fx=((x1-2.*x2+x3)*newcsilog+(x3-x1)/2.)
                    fy=((y1-2.*y2+y3)*newcsilog+(y3-y1)/2.)
                    jac=sqrt((fx**2)+(fy**2));                              !Jacobian           
                    intsg1A=intsg1A+((log((1.-csino)/2.)-1.)*NF(1)*jac*Jcsi*weightslog(m))
                    intsg2A=intsg2A+((log((1.-csino)/2.)-1.)*NF(2)*jac*Jcsi*weightslog(m))
                    intsg3A=intsg3A+((log((1.-csino)/2.)-1.)*NF(3)*jac*Jcsi*weightslog(m))
                END DO
                intnsg1B=0; intnsg2B=0; intnsg3B=0;  !Regular Gauss integration B - Before xd
                DO m=1,npg !-1 to 1
                    Jcsi=(csino+1.)/2.                                !Qsi jacobian
                    newcsi=Jcsi*csi(m)+(csino-1.)/2.                  !Qsi(eta)
                    k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2);
                    NF(1)=(newcsi-k2)*newcsi*k2/k3                    !Shape function of functional node 1 
                    NF(2)=(k1+k2)*(k1*k2+(k2-k1)*newcsi-newcsi**2)/k3 !Shape function of functional node 2
                    NF(3)=(k1+newcsi)*newcsi*k1/k3                    !Shape function of functional node 3
                    dxdcsi=((x1-2.*x2+x3)*newcsi+(x3-x1)/2.)
                    dydcsi=((y1-2.*y2+y3)*newcsi+(y3-y1)/2.)
                    jac=sqrt((dxdcsi**2)+(dydcsi**2));                !Jacobian      
                    RB1=((x3-x1)+(newcsi+csino)*(x1-2.*x2+x3))**2
                    RB2=((y3-y1)+(newcsi+csino)*(y1-2.*y2+y3))**2
                    RB=sqrt(RB1+RB2);                                 !RB
                    intnsg1B=intnsg1B+(log(RB)*NF(1)*jac*Jcsi*weights(m));
                    intnsg2B=intnsg2B+(log(RB)*NF(2)*jac*Jcsi*weights(m));
                    intnsg3B=intnsg3B+(log(RB)*NF(3)*jac*Jcsi*weights(m));
                END DO
                intsg1B=0; intsg2B=0; intsg3B=0; !Logarithmic Gauss integration B - Before xd
                DO m=1,npglog !0 to 1
                    Jcsi=(csino+1.)                                         !Qsi jacobian
                    newcsilog=-Jcsi*csilog(m)+csino                         !Qsi(eta)
                    k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2);
                    NF(1)=(newcsilog-k2)*newcsilog*k2/k3                    !Shape function of functional node 1
                    NF(2)=(k1+k2)*(k1*k2+(k2-k1)*newcsilog-newcsilog**2)/k3 !Shape function of functional node 2
                    NF(3)=(k1+newcsilog)*newcsilog*k1/k3                    !Shape function of functional node 3
                    fx=(((x1-2.*x2+x3)*newcsilog)+((x3-x1)/2.)) 
                    fy=(((y1-2.*y2+y3)*newcsilog)+((y3-y1)/2.)) 
                    jac=((fx**2)+(fy**2))**0.5;                             !Jacobian
                    intsg1B=intsg1B+((log((1.+csino)/2.)-1.)*NF(1)*jac*Jcsi*weightslog(m))
                    intsg2B=intsg2B+((log((1.+csino)/2.)-1.)*NF(2)*jac*Jcsi*weightslog(m))
                    intsg3B=intsg3B+((log((1.+csino)/2.)-1.)*NF(3)*jac*Jcsi*weightslog(m))
                END DO
                G_temp(i,s)=-(intnsg1A+intsg1A+intnsg1B+intsg1B)/(2*pi*kmaterial)
                G_temp(i,s+1)=-(intnsg2A+intsg2A+intnsg2B+intsg2B)/(2*pi*kmaterial)
                G_temp(i,s+2)=-(intnsg3A+intsg3A+intnsg3B+intsg3B)/(2*pi*kmaterial)
                H_temp(i,s)=0;
                H_temp(i,s+1)=0;
                H_temp(i,s+2)=0;                
            END IF
            s=s+3;
        END DO
    END DO
    DO i=1,ne*3 !Building H(i,i)
        H_temp(i,i)=0;
        DO j=1,ne*3
            IF (i.NE.j) THEN
                H_temp(i,i)=H_temp(i,i)-H_temp(i,j)
            END IF    
        END DO
    END DO 
    PRINT*," -Done."
    PRINT*, "  "

    !Building the known vector and auxiliary vector to build [H] and [G] 
    DO i=1,num_arestas_totais
        j=aresta_cc_temp(i)
        IF (j.NE.0)THEN
            termoconhearesta(j)=valor_temp_aresta(i)
            inversorvector1(j)=aresta_cc_temp(i)
        ELSE
            termoconhearesta(j)=0
            inversorvector1(j)=0
        END IF
    END DO
    DO i=1,num_arestas_totais
        j=aresta_cc_fluxo(i)
        IF (j.NE.0)THEN
            termoconhearesta(j)=valor_fluxo_aresta(i)
        ELSE
            termoconhearesta(j)=termoconhearesta(j)
        END IF
    END DO
    k=0;
    DO i=1,num_arestas_totais 
        j=num_ele_aresta(i)*3
        DO m=1,j 
            k=k+1
            termoconhe(k)=termoconhearesta(i)
            inversorvector2(k)=inversorvector1(i) !It's to store the known temperature node 
        END DO
    END DO
    DO i=1,ne*3 
        !Building [H] and [G] to use in solving the system of equation
        DO j=1,ne*3 !To variate the columns of [H] or [G]
            IF (inversorvector2(j).NE.0) THEN
                H(i,j)=(-1.)*G_temp(i,j)
                G(i,j)=(-1.)*H_temp(i,j)
            ELSE
                H(i,j)=H_temp(i,j)
                G(i,j)=G_temp(i,j)
            END IF
        END DO
        !Building the {b} vector to use in solving the system of equation
        btemp=0; 
        DO j=1,ne*3
            btemp=btemp+(G(i,j)*termoconhe(j))
        END DO
        b(i)=btemp
    END DO

    !QR decomposition method to solve the system of equation
    A(:,:)=H(:,:);
    PRINT*,"QR decomposition method to solve the system of equation."
    PRINT*," -Solving."
    DO k=1,ne*3,1                                           !To read the column
        IF (k.EQ.1) THEN                                    !To make the first column
            coluna_atual(:)=A(:,k)
            somanorma=0
            DO i=1,ne*3,1
                somanorma=somanorma+(coluna_atual(i)**2)
            END DO
            normau=sqrt(somanorma)                          !To normalize {u}
            Q(:,k)=coluna_atual(:)/normau
        ELSE
            proatual(:)=0
            coluna_atual(:)=A(:,k)                          !To make the others columns
            DO i=1,k-1,1                                    !Projection sum
                prodinter=sum(coluna_atual(:)*Q(:,i))       !Sum of dot product
                proatual(:)=proatual(:)+(prodinter*Q(:,i))
            END DO
            u(:)=coluna_atual(:)-proatual(:)                !Projection subtraction
            somanorma=sum(u(:)*u(:));
            normau=sqrt(somanorma)
            Q(:,k)=u(:)/normau                              !To write [Q]
        END IF
    END DO
    QT(:,:)=transpose(Q(:,:))                               !To transpose [Q]
    Rsystem=matmul(QT,A)                                    !To write [R]
    y=matmul(QT,b)
    x(ne*3)=y(ne*3)/Rsystem(ne*3,ne*3)
    DO i=(ne*3)-1,1,-1
        somats=0
        DO j=i+1,ne*3,1
            somats=somats+Rsystem(i,j)*x(j)
        END DO
        x(i)=(y(i)-somats)/Rsystem(i,i)
    END DO
    PRINT*, " -Done."
    PRINT*, "  "

    !Building {T} and {q} vectors
    DO i=1,ne*3
        IF (inversorvector2(i).EQ.0) THEN 
            Tbound(i)=x(i)          !{T} on boundary 
            Qbound(i)=termoconhe(i) !{q} on boundary
        ELSE
            Tbound(i)=termoconhe(i)
            Qbound(i)=x(i) 
        END IF
    END DO

    !Solving the internal points
    Hd(i,j)=0; Hdx(i,j)=0; Hdy(i,j)=0; Gd(i,j)=0; Gdx(i,j)=0; Gdy(i,j)=0
    DO i=1,neinternal !To define the source point 
        xd=xinternal(i)
        yd=yinternal(i)
        s=1;
        DO j=1,ne !To define the element
            x1=xelefinal(s); x2=xelefinal(s+1); x3=xelefinal(s+2)
            y1=yelefinal(s); y2=yelefinal(s+1); y3=yelefinal(s+2)
            intG1=0; intG2=0; intG3=0; intH1=0; intH2=0; intH3=0; 
            DO m=1,npg !Gauss regular integration to internal points
                FF1=(csi(m)**2-csi(m))/2.      
                FF2=1.-csi(m)**2               
                FF3=(csi(m)**2+csi(m))/2.      
                k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2); !Node shape function 
                NF(1)=(csi(m)-k2)*csi(m)*k2/k3
                NF(2)=(k1+k2)*(k1*k2+(k2-k1)*csi(m)-csi(m)**2)/k3      
                NF(3)=(k1+csi(m))*csi(m)*k1/k3  
                dFF1=(2.*csi(m)-1)/2.  
                dFF2=-2.*csi(m)         
                dFF3=(2.*csi(m)+1)/2.  
                xpi=(FF1*x1)+(FF2*x2)+(FF3*x3);            
                ypi=(FF1*y1)+(FF2*y2)+(FF3*y3);             
                Rx=(xpi-xd)
                Ry=(ypi-yd)
                R=((Rx**2)+(Ry**2))**0.5              !Distance between source point and field point (Radius) 
                fx=((x1-2.*x2+x3)*csi(m)+(x3-x1)/2.); 
                fy=((y1-2.*y2+y3)*csi(m)+(y3-y1)/2.);
                jac=((fx**2)+(fy**2))**0.5;           !Jacobian
                nx=fy/jac;                            !X component of the normal vector 
                ny=-fx/jac;                           !Y component of the normal vector 
                Tsolfun=(-1./(2.*pi*kmaterial))*log(R);                    !Temperature solfun
                Qsolfun=(1./(2.*pi*kmaterial))*((Rx*nx)+(Ry*ny))/(R**2);   !Flux solfun
                intG1=intG1+((Tsolfun*(NF(1))*jac*weights(m)));
                intG2=intG2+((Tsolfun*(NF(2))*jac*weights(m)));
                intG3=intG3+((Tsolfun*(NF(3))*jac*weights(m)));
                intH1=intH1+((Qsolfun*(NF(1))*jac*weights(m)));
                intH2=intH2+((Qsolfun*(NF(2))*jac*weights(m)));
                intH3=intH3+((Qsolfun*(NF(3))*jac*weights(m)));
            END DO
            Hd(i,s)=intH1
            Hd(i,s+1)=intH2
            Hd(i,s+2)=intH3
            Gd(i,s)=intG1
            Gd(i,s+1)=intG2
            Gd(i,s+2)=intG3
            !!!!!!
            intG1=0; intG2=0; intG3=0; intH1=0; intH2=0; intH3=0;
            DO m=1,npg !Gauss regular integration
                FF1=(csi(m)**2-csi(m))/2.      
                FF2=1.-csi(m)**2               
                FF3=(csi(m)**2+csi(m))/2.      
                k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2); !Node shape function 
                NF(1)=(csi(m)-k2)*csi(m)*k2/k3
                NF(2)=(k1+k2)*(k1*k2+(k2-k1)*csi(m)-csi(m)**2)/k3      
                NF(3)=(k1+csi(m))*csi(m)*k1/k3  
                dFF1=(2.*csi(m)-1)/2.  
                dFF2=-2.*csi(m)         
                dFF3=(2.*csi(m)+1)/2.  
                xpi=(FF1*x1)+(FF2*x2)+(FF3*x3);            
                ypi=(FF1*y1)+(FF2*y2)+(FF3*y3);             
                Rx=(xpi-xd)
                Ry=(ypi-yd)
                R=((Rx**2)+(Ry**2))**0.5              !Distance between source point and field point (Radius) 
                fx=((x1-2.*x2+x3)*csi(m)+(x3-x1)/2.); 
                fy=((y1-2.*y2+y3)*csi(m)+(y3-y1)/2.);
                jac=((fx**2)+(fy**2))**0.5;           !Jacobian
                dx=x1*dFF1+x2*dFF2+x3*dFF3
                dy=y1*dFF1+y2*dFF2+y3*dFF3
                nx=dy/jac;                            !X component of the normal vector 
                ny=-dx/jac;                           !Y component of the normal vector 
                intG1=intG1+(((1./(2.*pi*kmaterial))*(NF(1))*(Rx/(R**2))*jac*weights(m)));
                intG2=intG2+(((1./(2.*pi*kmaterial))*(NF(2))*(Rx/(R**2))*jac*weights(m)));
                intG3=intG3+(((1./(2.*pi*kmaterial))*(NF(3))*(Rx/(R**2))*jac*weights(m)));
                intH1=intH1+((-1./(2.*pi*kmaterial))*(NF(1))*(2*Rx*Ry*ny+(Rx**2-Ry**2)*nx)*&
                      jac*weights(m)/(R**4));
                intH2=intH2+((-1./(2.*pi*kmaterial))*(NF(2))*(2*Rx*Ry*ny+(Rx**2-Ry**2)*nx)*&
                      jac*weights(m)/(R**4));
                intH3=intH3+((-1./(2.*pi*kmaterial))*(NF(3))*(2*Rx*Ry*ny+(Rx**2-Ry**2)*nx)*&
                      jac*weights(m)/(R**4));
            END DO
            Hdx(i,s)=intH1
            Hdx(i,s+1)=intH2
            Hdx(i,s+2)=intH3
            Gdx(i,s)=intG1
            Gdx(i,s+1)=intG2
            Gdx(i,s+2)=intG3
            !!!!!!
            intG1=0; intG2=0; intG3=0; intH1=0; intH2=0; intH3=0;
            DO m=1,npg !Gauss regular integration
                FF1=(csi(m)**2-csi(m))/2.      
                FF2=1.-csi(m)**2               
                FF3=(csi(m)**2+csi(m))/2.      
                k1=2./3.; k2=2./3.; k3=k1*k2*(k1+k2); !Node shape function 
                NF(1)=(csi(m)-k2)*csi(m)*k2/k3
                NF(2)=(k1+k2)*(k1*k2+(k2-k1)*csi(m)-csi(m)**2)/k3      
                NF(3)=(k1+csi(m))*csi(m)*k1/k3  
                dFF1=(2.*csi(m)-1)/2.  
                dFF2=-2.*csi(m)         
                dFF3=(2.*csi(m)+1)/2.  
                xpi=(FF1*x1)+(FF2*x2)+(FF3*x3);            
                ypi=(FF1*y1)+(FF2*y2)+(FF3*y3);             
                Rx=(xpi-xd)
                Ry=(ypi-yd)
                R=((Rx**2)+(Ry**2))**0.5              !Distance between source point and field point (Radius) 
                fx=((x1-2.*x2+x3)*csi(m)+(x3-x1)/2.); 
                fy=((y1-2.*y2+y3)*csi(m)+(y3-y1)/2.);
                jac=((fx**2)+(fy**2))**0.5;           !Jacobian
                dx=x1*dFF1+x2*dFF2+x3*dFF3
                dy=y1*dFF1+y2*dFF2+y3*dFF3
                nx=dy/jac;                            !X component of the normal vector 
                ny=-dx/jac;                           !Y component of the normal vector 
                intG1=intG1+(((1./(2.*pi*kmaterial))*(NF(1))*(Ry/(R**2))*jac*weights(m)));
                intG2=intG2+(((1./(2.*pi*kmaterial))*(NF(2))*(Ry/(R**2))*jac*weights(m)));
                intG3=intG3+(((1./(2.*pi*kmaterial))*(NF(3))*(Ry/(R**2))*jac*weights(m)));
                intH1=intH1+((-1./(2.*pi*kmaterial))*(NF(1))*(2*Rx*Ry*nx+(Ry**2-Rx**2)*ny)*&
                      jac*weights(m)/(R**4));
                intH2=intH2+((-1./(2.*pi*kmaterial))*(NF(2))*(2*Rx*Ry*nx+(Ry**2-Rx**2)*ny)*&
                      jac*weights(m)/(R**4));
                intH3=intH3+((-1./(2.*pi*kmaterial))*(NF(3))*(2*Rx*Ry*nx+(Ry**2-Rx**2)*ny)*&
                      jac*weights(m)/(R**4));
            END DO
            Hdy(i,s)=intH1
            Hdy(i,s+1)=intH2
            Hdy(i,s+2)=intH3
            Gdy(i,s)=intG1
            Gdy(i,s+1)=intG2
            Gdy(i,s+2)=intG3
            !!!!!!
            s=s+3;
        END DO
    END DO
    DO i=1,neinternal
        Tin(i)=0; Qinx(i)=0; Qiny(i)=0;
        DO j=1,ne*3
            Tin(i)=Tin(i)-Gd(i,j)*Qbound(j)+Hd(i,j)*Tbound(j)
            Qinx(i)=Qinx(i)+Gdx(i,j)*Qbound(j)+Hdx(i,j)*Tbound(j)
            Qiny(i)=Qiny(i)+Gdy(i,j)*Qbound(j)+Hdy(i,j)*Tbound(j)
        END DO
    END DO
   
    !OUTPUT DATA
    PRINT*,"VTK file writer."
    PRINT*," -Writing."

    OPEN(unit=200,file='resultsq.vtk',status='old') !Ploting intern flux
    WRITE(200,*)"x,y,xvector,yvector"
    DO i=1,neinternal 
        WRITE(200,*)xinternal(i),",",yinternal(i),",",Qinx(i),",",Qiny(i)
    END DO
    OPEN(unit=300,file='resultst.vtk',status='old') !Ploting temperature
    DO i=1,ne*3
        WRITE(300,*)xnofinal(i),",",ynofinal(i),",",0.,",",Tbound(i)
    END DO
    DO i=1,neinternal
        WRITE(300,*)xinternal(i),",",yinternal(i),",",0.,",",Tin(i)
    END DO
    OPEN(unit=400,file='resultsqbound.vtk',status='old') !Ploting flux on boundary
    DO i=1,ne*3
        WRITE(400,*)xnofinal(i),",",ynofinal(i),",",0.,",",Qbound(i)
    END DO
    PRINT*," -Done."
    
    !Ending timer
    CALL system_clock(count_1, count_rate, count_max)
    time_final=count_1*1.0/count_rate
    !Elapsed timer
    elapsed_time=time_final-time_init
    ! Write elapsed time
    PRINT*, " "
    PRINT*, "Elapsed time:",elapsed_time-int(elapsed_time),"seconds"
    PRINT*, "Developed by Caio Moura"
    PRINT*, "FEM-UNICAMP, 2021"

END PROGRAM