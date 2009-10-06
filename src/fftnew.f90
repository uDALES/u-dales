!> \file fftnew.f90
!! Performs Fast Fourier Transforms.
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

  subroutine RADB2 (IDO,L1,CC,CH,WA1)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2), WA1(IDO-1)
  do K=1,L1
     CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
     CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  end do

  if (ido > 2 ) then
    IDP2 = IDO+2
    do K=1,L1
      do I=3,IDO,2
        IC = IDP2-I
        CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
        TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
        CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
        TI2 = CC(I,1,K)+CC(IC,2,K)
        CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
        CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
      end do
    end do
    if (MOD(IDO,2) == 1) return
  end if
  if (ido >=2) then
    do K=1,L1
       CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
       CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
    end do
  end if
  return
  end

  subroutine RADB3 (IDO,L1,CC,CH,WA1,WA2)
      DIMENSION       CC(IDO,3,L1) ,CH(IDO,L1,3) ,WA1(1)     ,WA2(1)
  DATA TAUR,TAUI /-.5,.866025403784439/
  do K=1,L1
     TR2 = CC(IDO,2,K)+CC(IDO,2,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
  end do
  if (IDO == 1) return
  IDP2 = IDO+2
  do K=1,L1
     do I=3,IDO,2
        IC = IDP2-I
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,3,K)-CC(IC,2,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
        CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
      end do
    end do
  return
  end

  subroutine RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      DIMENSION       CC(IDO,4,L1),CH(IDO,L1,4),WA1(IDO-1),WA2(IDO-1),WA3(IDO-1)
  DATA SQRT2 /1.414213562373095/
  do K=1,L1
     TR1 = CC(1,1,K)-CC(IDO,4,K)
     TR2 = CC(1,1,K)+CC(IDO,4,K)
     TR3 = CC(IDO,2,K)+CC(IDO,2,K)
     TR4 = CC(1,3,K)+CC(1,3,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,2) = TR1-TR4
     CH(1,K,3) = TR2-TR3
     CH(1,K,4) = TR1+TR4
  end do
  if (ido>2) then
    IDP2 = IDO+2
    do K=1,L1
       do I=3,IDO,2
          IC = IDP2-I
          TI1 = CC(I,1,K)+CC(IC,4,K)
          TI2 = CC(I,1,K)-CC(IC,4,K)
          TI3 = CC(I,3,K)-CC(IC,2,K)
          TR4 = CC(I,3,K)+CC(IC,2,K)
          TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
          TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
          TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
          TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
          CH(I-1,K,1) = TR2+TR3
          CR3 = TR2-TR3
          CH(I,K,1) = TI2+TI3
          CI3 = TI2-TI3
          CR2 = TR1-TR4
          CR4 = TR1+TR4
          CI2 = TI1+TI4
          CI4 = TI1-TI4
          CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
          CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
          CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
          CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
          CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
          CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
        end do
    end do
    if (MOD(IDO,2) == 1) return
  endif
  if (ido>=2) then
    do K=1,L1
       TI1 = CC(1,2,K)+CC(1,4,K)
       TI2 = CC(1,4,K)-CC(1,2,K)
       TR1 = CC(IDO,1,K)-CC(IDO,3,K)
       TR2 = CC(IDO,1,K)+CC(IDO,3,K)
       CH(IDO,K,1) = TR2+TR2
       CH(IDO,K,2) = SQRT2*(TR1-TI1)
       CH(IDO,K,3) = TI2+TI2
       CH(IDO,K,4) = -SQRT2*(TR1+TI1)
    end do
  endif
  return
  end

  subroutine RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DIMENSION  CC(IDO,5,L1),CH(IDO,L1,5),WA1(1),WA2(1),WA3(1),WA4(1)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,-.809016994374947,.587785252292473/
  do K=1,L1
     TI5 = CC(1,3,K)+CC(1,3,K)
     TI4 = CC(1,5,K)+CC(1,5,K)
     TR2 = CC(IDO,2,K)+CC(IDO,2,K)
     TR3 = CC(IDO,4,K)+CC(IDO,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI5 = TI11*TI5+TI12*TI4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(1,K,5) = CR2+CI5
  end do
  if (IDO == 1) return
  IDP2 = IDO+2
  do K=1,L1
    do I=3,IDO,2
        IC = IDP2-I
        TI5 = CC(I,3,K)+CC(IC,2,K)
        TI2 = CC(I,3,K)-CC(IC,2,K)
        TI4 = CC(I,5,K)+CC(IC,4,K)
        TI3 = CC(I,5,K)-CC(IC,4,K)
        TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
        TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
        CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
        CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
        CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
        CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
     end do
  end do
  return
  end

  subroutine RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DIMENSION  CH(IDO,L1,IP),CC(IDO,IP,L1),C1(IDO,L1,IP),C2(IDL1,IP), CH2(IDL1,IP) ,WA(1)
  DATA TPI/6.28318530717959/
  ARG = TPI/FLOAT(IP)
  DCP = COS(ARG)
  DSP = SIN(ARG)
  IDP2 = IDO+2
  NBD = (IDO-1)/2
  IPP2 = IP+2
  IPPH = (IP+1)/2
  if (IDO >= L1) then
    do K=1,L1
    do I=1,IDO
      CH(I,K,1) = CC(I,1,K)
    end do
    end do
  else
    do I=1,IDO
    do K=1,L1
      CH(I,K,1) = CC(I,1,K)
    end do
    end do
  endif
  do J=2,IPPH
     JC = IPP2-J
     J2 = J+J
    do K=1,L1
      CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
      CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
    end do
  end do
  if (IDO /= 1) then
    if (NBD >= L1) then
      do J=2,IPPH
        JC = IPP2-J
         do K=1,L1
            do I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
            end do
         end do
      end do
    else
      do J=2,IPPH
        JC = IPP2-J
        do I=3,IDO,2
          IC = IDP2-I
          do K=1,L1
            CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
            CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
            CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
            CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
          end do
        end do
      end do
    end if
  end if
  AR1 = 1.
  AI1 = 0.
  do L=2,IPPH
     LC = IPP2-L
     AR1H = DCP*AR1-DSP*AI1
     AI1 = DCP*AI1+DSP*AR1
     AR1 = AR1H
     do IK=1,IDL1
        C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
        C2(IK,LC) = AI1*CH2(IK,IP)
     end do
     DC2 = AR1
     DS2 = AI1
     AR2 = AR1
     AI2 = AI1
     do J=3,IPPH
        JC = IPP2-J
        AR2H = DC2*AR2-DS2*AI2
        AI2 = DC2*AI2+DS2*AR2
        AR2 = AR2H
        do IK=1,IDL1
           C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
           C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
        end do
      end do
  end do
  do J=2,IPPH
     do IK=1,IDL1
        CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
     end do
  end do
  do J=2,IPPH
     JC = IPP2-J
     do K=1,L1
        CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
        CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
     end do
   end do
  if (IDO /= 1) then
    if (NBD >= L1) then
      do J=2,IPPH
        JC = IPP2-J
        do K=1,L1
          do I=3,IDO,2
            CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
            CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
            CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
            CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
          end do
        end do
      end do
    else
      do J=2,IPPH
        JC = IPP2-J
        do I=3,IDO,2
          do K=1,L1
            CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
            CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
            CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
            CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
          end do
        end do
      end do
    end if
  end if
  if (IDO == 1) return
  do IK=1,IDL1
     C2(IK,1) = CH2(IK,1)
  end do
  do J=2,IP
     do K=1,L1
        C1(1,K,J) = CH(1,K,J)
     end do
  end do
  if (NBD <= L1) then
    IS = -IDO
    do J=2,IP
      IS = IS+IDO
      IDIJ = IS
      do I=3,IDO,2
        IDIJ = IDIJ+2
        do K=1,L1
          C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
          C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        end do
      end do
    end do
  else
    IS = -IDO
    do J=2,IP
      IS = IS+IDO
      do K=1,L1
        IDIJ = IS
        do I=3,IDO,2
          IDIJ = IDIJ+2
          C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
          C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
        end do
      end do
    end do
  end if
  return
  end

  subroutine RADF2 (IDO,L1,CC,CH,WA1)
      DIMENSION       CH(IDO,2,L1),CC(IDO,L1,2),WA1(IDO-1)
  do K=1,L1
     CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
     CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  end do
  if (ido >2) then
    IDP2 = IDO+2
    do K=1,L1
      do I=3,IDO,2
        IC = IDP2-I
        TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CH(I,1,K) = CC(I,K,1)+TI2
        CH(IC,2,K) = TI2-CC(I,K,1)
        CH(I-1,1,K) = CC(I-1,K,1)+TR2
        CH(IC-1,2,K) = CC(I-1,K,1)-TR2
      end do
    end do
    if (MOD(IDO,2) == 1) return
  end if
  if (ido >=2) then
    do K=1,L1
      CH(1,2,K) = -CC(IDO,K,2)
      CH(IDO,1,K) = CC(IDO,K,1)
    end do
  end if
  return
  end

  subroutine RADF3 (IDO,L1,CC,CH,WA1,WA2)
      DIMENSION       CH(IDO,3,L1),CC(IDO,L1,3), WA1(1)     ,WA2(1)
  DATA TAUR,TAUI /-.5,.866025403784439/
  do K=1,L1
     CR2 = CC(1,K,2)+CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2
     CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
     CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  end do
  if (IDO == 1) return
  IDP2 = IDO+2
  do K=1,L1
     do I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR2 = DR2+DR3
        CI2 = DI2+DI3
        CH(I-1,1,K) = CC(I-1,K,1)+CR2
        CH(I,1,K) = CC(I,K,1)+CI2
        TR2 = CC(I-1,K,1)+TAUR*CR2
        TI2 = CC(I,K,1)+TAUR*CI2
        TR3 = TAUI*(DI2-DI3)
        TI3 = TAUI*(DR3-DR2)
        CH(I-1,3,K) = TR2+TR3
        CH(IC-1,2,K) = TR2-TR3
        CH(I,3,K) = TI2+TI3
        CH(IC,2,K) = TI3-TI2
     end do
  end do
  return
  end

  subroutine RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      DIMENSION       CC(IDO,L1,4),CH(IDO,4,L1), WA1(IDO-1)     ,WA2(IDO-1)     ,WA3(IDO-1)
  DATA HSQT2 /.7071067811865475/
  do K=1,L1
     TR1 = CC(1,K,2)+CC(1,K,4)
     TR2 = CC(1,K,1)+CC(1,K,3)
     CH(1,1,K) = TR1+TR2
     CH(IDO,4,K) = TR2-TR1
     CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
     CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  end do
  if (ido>2) then
    IDP2 = IDO+2
    do K=1,L1
      do I=3,IDO,2
        IC = IDP2-I
        CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        TR1 = CR2+CR4
        TR4 = CR4-CR2
        TI1 = CI2+CI4
        TI4 = CI2-CI4
        TI2 = CC(I,K,1)+CI3
        TI3 = CC(I,K,1)-CI3
        TR2 = CC(I-1,K,1)+CR3
        TR3 = CC(I-1,K,1)-CR3
        CH(I-1,1,K) = TR1+TR2
        CH(IC-1,4,K) = TR2-TR1
        CH(I,1,K) = TI1+TI2
        CH(IC,4,K) = TI1-TI2
        CH(I-1,3,K) = TI4+TR3
        CH(IC-1,2,K) = TR3-TI4
        CH(I,3,K) = TR4+TI3
        CH(IC,2,K) = TR4-TI3
      end do
    end do
    if (MOD(IDO,2) == 1) return
  end if
  if (ido>=2) then
    do K=1,L1
      TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
      TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
      CH(IDO,1,K) = TR1+CC(IDO,K,1)
      CH(IDO,3,K) = CC(IDO,K,1)-TR1
      CH(1,2,K) = TI1-CC(IDO,K,3)
      CH(1,4,K) = TI1+CC(IDO,K,3)
    end do
  end if
  return
  end

  subroutine RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DIMENSION       CC(IDO,L1,5),CH(IDO,5,L1),WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,-.809016994374947,.587785252292473/
  do K=1,L1
     CR2 = CC(1,K,5)+CC(1,K,2)
     CI5 = CC(1,K,5)-CC(1,K,2)
     CR3 = CC(1,K,4)+CC(1,K,3)
     CI4 = CC(1,K,4)-CC(1,K,3)
     CH(1,1,K) = CC(1,K,1)+CR2+CR3
     CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
     CH(1,3,K) = TI11*CI5+TI12*CI4
     CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
     CH(1,5,K) = TI12*CI5-TI11*CI4
  end do
  if (IDO == 1) return
  IDP2 = IDO+2
  do K=1,L1
     do I=3,IDO,2
        IC = IDP2-I
        DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
        DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
        DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
        DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
        DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
        DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
        DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
        DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
        CR2 = DR2+DR5
        CI5 = DR5-DR2
        CR5 = DI2-DI5
        CI2 = DI2+DI5
        CR3 = DR3+DR4
        CI4 = DR4-DR3
        CR4 = DI3-DI4
        CI3 = DI3+DI4
        CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
        CH(I,1,K) = CC(I,K,1)+CI2+CI3
        TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
        TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
        TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
        TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
        TR5 = TI11*CR5+TI12*CR4
        TI5 = TI11*CI5+TI12*CI4
        TR4 = TI12*CR5-TI11*CR4
        TI4 = TI12*CI5-TI11*CI4
        CH(I-1,3,K) = TR2+TR5
        CH(IC-1,2,K) = TR2-TR5
        CH(I,3,K) = TI2+TI5
        CH(IC,2,K) = TI5-TI2
        CH(I-1,5,K) = TR3+TR4
        CH(IC-1,4,K) = TR3-TR4
        CH(I,5,K) = TI3+TI4
        CH(IC,4,K) = TI4-TI3
    end do
  end do
  return
  end

  subroutine RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DIMENSION CH(IDO,L1,IP),CC(IDO,IP,L1),C1(IDO,L1,IP),C2(IDL1,IP),CH2(IDL1,IP),WA(1)
  DATA TPI/6.28318530717959/
  ARG = TPI/FLOAT(IP)
  DCP = COS(ARG)
  DSP = SIN(ARG)
  IPPH = (IP+1)/2
  IPP2 = IP+2
  IDP2 = IDO+2
  NBD = (IDO-1)/2
  if (IDO /= 1) then
    do IK=1,IDL1
       CH2(IK,1) = C2(IK,1)
    end do
    do J=2,IP
      do K=1,L1
        CH(1,K,J) = C1(1,K,J)
      end do
    end do
    if (NBD <= L1) then
      IS = -IDO
      do J=2,IP
        IS = IS+IDO
        IDIJ = IS
        do I=3,IDO,2
          IDIJ = IDIJ+2
          do K=1,L1
            CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
            CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
          end do
        end do
      end do
    else
      IS = -IDO
      do J=2,IP
        IS = IS+IDO
        do K=1,L1
          IDIJ = IS
          do I=3,IDO,2
            IDIJ = IDIJ+2
            CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
            CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
          end do
        end do
      end do
    end if
    if (NBD >= L1) then
      do J=2,IPPH
        JC = IPP2-J
        do K=1,L1
          do I=3,IDO,2
            C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
            C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
            C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
            C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
          end do
        end do
      end do
    else
      do J=2,IPPH
        JC = IPP2-J
        do I=3,IDO,2
          do K=1,L1
            C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
            C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
            C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
            C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
          end do
        end do
      end do
    end if
  else
    do IK=1,IDL1
      C2(IK,1) = CH2(IK,1)
    end do
  end if
  do J=2,IPPH
    JC = IPP2-J
    do K=1,L1
      C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
      C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
    end do
  end do

  AR1 = 1.
  AI1 = 0.
  do L=2,IPPH
    LC = IPP2-L
    AR1H = DCP*AR1-DSP*AI1
    AI1 = DCP*AI1+DSP*AR1
    AR1 = AR1H
    do IK=1,IDL1
      CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
      CH2(IK,LC) = AI1*C2(IK,IP)
    end do
    DC2 = AR1
    DS2 = AI1
    AR2 = AR1
    AI2 = AI1
    do J=3,IPPH
      JC = IPP2-J
      AR2H = DC2*AR2-DS2*AI2
      AI2 = DC2*AI2+DS2*AR2
      AR2 = AR2H
      do IK=1,IDL1
        CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
        CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
      end do
    end do
  end do
  do J=2,IPPH
    do IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+C2(IK,J)
    end do
  end do

  if (IDO >= L1) then
    do K=1,L1
      do I=1,IDO
        CC(I,1,K) = CH(I,K,1)
      end do
    end do
  else
    do I=1,IDO
      do K=1,L1
        CC(I,1,K) = CH(I,K,1)
      end do
    end do
  end if
  do J=2,IPPH
    JC = IPP2-J
    J2 = J+J
    do K=1,L1
      CC(IDO,J2-2,K) = CH(1,K,J)
      CC(1,J2-1,K) = CH(1,K,JC)
    end do
  end do
  if (IDO /= 1) then
  if (NBD >= L1) then
    do J=2,IPPH
      JC = IPP2-J
      J2 = J+J
      do K=1,L1
        do I=3,IDO,2
          IC = IDP2-I
          CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
          CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
          CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
          CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
        end do
      end do
    end do
  else
    do J=2,IPPH
      JC = IPP2-J
      J2 = J+J
      do I=3,IDO,2
        IC = IDP2-I
        do K=1,L1
          CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
          CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
          CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
          CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
        end do
      end do
    end do
  end if
  end if
  return
  end

  subroutine RFFTB (N,R,WSAVE)
      DIMENSION       R(1)       ,WSAVE(2*N+1)
  if (N == 1) return
  call RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
  return
  end

  subroutine RFFTB1 (N,C,CH,WA,IFAC)
      DIMENSION       CH(N)      ,C(N)       ,WA(N+1)      ,IFAC(2*N+1)
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  do K1=1,NF
    IP = IFAC(K1+2)
    L2 = IP*L1
    IDO = N/L2
    IDL1 = IDO*L1

    select case (IP)

    case (4)
      IX2 = IW+IDO
      IX3 = IX2+IDO
      if (NA == 0) then
        call RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      else
        call RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      end if
      NA = 1-NA

    case (2)
      if (NA == 0) then
        call RADB2 (IDO,L1,C,CH,WA(IW))
      else
        call RADB2 (IDO,L1,CH,C,WA(IW))
      end if
      NA = 1-NA

    case (3)
      IX2 = IW+IDO
      if (NA == 0) then
        call RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
      else
        call RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
      end if
      NA = 1-NA

    case (5)
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IX4 = IX3+IDO
      if (NA == 0) then
        call RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      else
        call RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      end if
      NA = 1-NA

    case default
      if (NA == 0) then
        call RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
      else
        call RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
      end if
      if (IDO == 1) NA = 1-NA
    end select

    L1 = L2
    IW = IW+(IP-1)*IDO
  end do
  if (NA == 0) return
  do I=1,N
     C(I) = CH(I)
  end do
  return
  end

  subroutine RFFTF (N,R,WSAVE)
      DIMENSION       R(1)       ,WSAVE(2*N+1)
  if (N == 1) return
  call RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
  return
  end

  subroutine RFFTF1 (N,C,CH,WA,IFAC)
      DIMENSION       CH(N)      ,C(N)       ,WA(N+1)      ,IFAC(2*N+1)
  NF = IFAC(2)
  NA = 1
  L2 = N
  IW = N
  do K1=1,NF
    KH = NF-K1
    IP = IFAC(KH+3)
    L1 = L2/IP
    IDO = N/L2
    IDL1 = IDO*L1
    IW = IW-(IP-1)*IDO
    NA = 1-NA

    select case (IP)

    case (4)
      IX2 = IW+IDO
      IX3 = IX2+IDO
      if (NA == 0) then
        call RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      else
        call RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      end if

    case (2)
      if (NA == 0) then
        call RADF2 (IDO,L1,C,CH,WA(IW))
      else
        call RADF2 (IDO,L1,CH,C,WA(IW))
      end if

    case(3)
      IX2 = IW+IDO
      if (NA == 0) then
        call RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
      else
        call RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
      end if

    case (5)
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IX4 = IX3+IDO
      if (NA == 0) then
        call RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      else
        call RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      end if

    case default
      if (IDO == 1) NA = 1-NA
      if (NA == 0) then
        call RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
        NA = 1
      else
        call RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
        NA = 0
      end if
    end select
    L2 = L1
  end do
  if (NA == 1) return
  do I=1,N
     C(I) = CH(I)
  end do
  return
  end

  subroutine RFFTI (N,WSAVE)
      DIMENSION       WSAVE(2*N+1)
  if (N == 1) return
  call RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
  return
  end

  subroutine RFFTI1 (N,WA,IFAC)
  logical firstloop
  DIMENSION       WA(N+1)      ,IFAC(2*N+1)    ,NTRYH(4)
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/

  NTRY = 4
  NQ   = 0
  NL = N
  NF = 0
  NR = 1
  J = 0
  firstloop = .true.
  do
    do while (NR /= 0)
      if (firstloop) then
        J = J+1
        if (J<=4) then
          NTRY = NTRYH(J)
        else
          NTRY = NTRY+2
        endif
      else
        firstloop = .true.
      endif
      NQ = NL/NTRY
      NR = NL-NTRY*NQ
    end do
    NF = NF+1
    IFAC(NF+2) = NTRY
    NL = NQ
    if ((NTRY == 2) .and. (NF /=1)) then
      do I=2,NF
        IB = NF-I+2
        IFAC(IB+2) = IFAC(IB+1)
      end do
      IFAC(3) = 2
    end if
    if (NL == 1) exit
    firstloop = .false.
    NR = 1
  end do
  IFAC(1) = N
  IFAC(2) = NF
  TPI = 6.28318530717959
  ARGH = TPI/FLOAT(N)
  IS = 0
  NFM1 = NF-1
  L1 = 1
  if (NFM1 == 0) return
  do K1=1,NFM1
    IP = IFAC(K1+2)
    LD = 0
    L2 = L1*IP
    IDO = N/L2
    IPM = IP-1
    do J=1,IPM
      LD = LD+L1
      I = IS
      ARGLD = FLOAT(LD)*ARGH
      FI = 0.
      do II=3,IDO,2
        I = I+2
        FI = FI+1.
        ARG = FI*ARGLD
        WA(I-1) = COS(ARG)
        WA(I) = SIN(ARG)
      end do
      IS = IS+IDO
    end do
    L1 = L2
  end do
  return
  end
