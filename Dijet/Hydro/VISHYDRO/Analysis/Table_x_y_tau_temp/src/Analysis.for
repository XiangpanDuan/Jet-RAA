      Program Main

        Implicit None
        Double Precision:: tau,x,y,e,s,Temp,Bd,vx,vy,Temp0
        Integer:: itau, ix, iy

        Open(66, FILE='../Output/x_y_tau_temp.dat', STATUS='UNKNOWN')

        Call hydroSetFiles()   ! set internal data file reader
        Call hydroReadInfo2D(0.6d0,0d0,0d0,e,s,Temp0,Bd,vx,vy)
        Write(*,*) Temp0;
        Do ix = -150, 150, 1
          Do iy = -150, 150, 1
            Do itau = 0, 200, 1
              tau=itau*0.1d0;
              x=ix*0.1d0
              y=iy*0.1d0
              Call hydroReadInfo2D(tau,x,y,e,s,Temp,Bd,vx,vy)
              Write(66,'(3(F8.4, 2x), 2x, F20.18)') x, y, tau, Temp
            EndDo
          EndDo
        EndDo


      End Program