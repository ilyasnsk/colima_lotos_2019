

SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
common/nkr_max/nmax
real U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
DO J=1,N
	S=0.
	IF(W(J).NE.0.)THEN
		DO I=1,M
			S=S+U(I,J)*B(I)
		end do
		S=S/W(J)
	ENDIF
	TMP(J)=S
end do
DO J=1,N
	S=0.
	DO JJ=1,N
		S=S+V(J,JJ)*TMP(JJ)
	end do
	X(J)=S
end do
RETURN
END


SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
common/nkr_max/nmax


real A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)

G=0.0
SCALE=0.0
ANORM=0.0
DO I=1,N
	L=I+1
	RV1(I)=SCALE*G
	G=0.0
	S=0.0
	SCALE=0.0

	IF (I.LE.M) THEN
		do K=I,M
			SCALE=SCALE+ABS(A(K,I))
		end do
		IF (SCALE.NE.0.0) THEN
			DO K=I,M
				A(K,I)=A(K,I)/SCALE
				S=S+A(K,I)*A(K,I)
			end do
			F=A(I,I)
			G=-SIGN(SQRT(S),F)
			H=F*G-S
			A(I,I)=F-G
			IF (I.NE.N) THEN
				DO J=L,N
					S=0.0
					DO K=I,M
						S=S+A(K,I)*A(K,J)
					end do
					F=S/H
					DO K=I,M
						A(K,J)=A(K,J)+F*A(K,I)
					end do
				end do
			ENDIF

			DO K= I,M
			  A(K,I)=SCALE*A(K,I)
			end do
		ENDIF
	ENDIF

	W(I)=SCALE *G
	G=0.0
	S=0.0
	SCALE=0.0
	IF ((I.LE.M).AND.(I.NE.N)) THEN
		DO K=L,N
			SCALE=SCALE+ABS(A(I,K))
		end do
		IF (SCALE.NE.0.0) THEN
			DO K=L,N
				A(I,K)=A(I,K)/SCALE
				S=S+A(I,K)*A(I,K)
			end do
			F=A(I,L)
			G=-SIGN(SQRT(S),F)
			H=F*G-S
			A(I,L)=F-G
			DO K=L,N
				RV1(K)=A(I,K)/H
			end do
			IF (I.NE.M) THEN
			DO J=L,M
				S=0.0
				DO K=L,N
					S=S+A(J,K)*A(I,K)
				end do
				DO K=L,N
					A(J,K)=A(J,K)+S*RV1(K)
				end do
			end do
			ENDIF
			DO K=L,N
			A(I,K)=SCALE*A(I,K)
			end do
		ENDIF
	ENDIF
	ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
end do

DO I=N,1,-1
	IF (I.LT.N) THEN
		IF (G.NE.0.0) THEN
			DO J=L,N
				V(J,I)=(A(I,J)/A(I,L))/G
			end do
			DO J=L,N
				S=0.0

				DO K=L,N
					S=S+A(I,K)*V(K,J)
				end do

				DO K=L,N
					V(K,J)=V(K,J)+S*V(K,I)
				end do

			end do
		ENDIF
		DO J=L,N
			V(I,J)=0.0
			V(J,I)=0.0
		end do
	ENDIF
	V(I,I)=1.0
	G=RV1(I)
	L=I
end do



DO I=N,1,-1
	L=I+1
	G=W(I)
	IF (I.LT.N) THEN
		DO J=L,N
			A(I,J)=0.0
		end do
	ENDIF
	IF (G.NE.0.0) THEN
		G=1.0/G
		IF (I.NE.N) THEN
			DO J=L,N
				S=0.0
				DO K=L,M
					S=S+A(K,I)*A(K,J)
				end do
				F=(S/A(I,I))*G
				DO K=I,M
					A(K,J)=A(K,J)+F*A(K,I)
				end do
			end do
		ENDIF
		DO J=I,M
			A(J,I)=A(J,I)*G
		end do
	ELSE
		DO J= I,M
			A(J,I)=0.0
		end do
	ENDIF
	A(I,I)=A(I,I)+1.0
end do



DO K=N,1,-1
	DO ITS=1,30
		DO L=K,1,-1
			NM=L-1
			IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
			IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
		end do


1       C=0.0
		S=1.0
		DO I=L,K
			F=S*RV1(I)
			IF ((ABS(F)+ANORM).NE.ANORM) THEN
				G=W(I)
				H=SQRT(F*F+G*G)
				W(I)=H
				H=1.0/H
				C= (G*H)
				S=-(F*H)
				DO J=1,M
					Y=A(J,NM)
					Z=A(J,I)
					A(J,NM)=(Y*C)+(Z*S)
					A(J,I)=-(Y*S)+(Z*C)
				end do
			ENDIF
		end do


2       Z=W(K)
		IF (L.EQ.K) THEN
			IF (Z.LT.0.0) THEN
				W(K)=-Z
				DO J=1,N
					V(J,K)=-V(J,K)
				end do
			ENDIF
			GO TO 3
		ENDIF
		!         IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
		X=W(L)
		NM=K-1
		Y=W(NM)
		G=RV1(NM)
		H=RV1(K)
		F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
		G=SQRT(F*F+1.0)
		F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
		C=1.0
		S=1.0
		DO J=L,NM
			I=J+1
			G=RV1(I)
			Y=W(I)
			H=S*G
			G=C*G
			Z=SQRT(F*F+H*H)
			RV1(J)=Z
			C=F/Z
			S=H/Z
			F= (X*C)+(G*S)
			G=-(X*S)+(G*C)
			H=Y*S
			Y=Y*C
			DO NM=1,N
				X=V(NM,J)
				Z=V(NM,I)
				V(NM,J)= (X*C)+(Z*S)
				V(NM,I)=-(X*S)+(Z*C)
			end do
			Z=SQRT(F*F+H*H)
			W(J)=Z
			IF (Z.NE.0.0) THEN
				Z=1.0/Z
				C=F*Z
				S=H*Z
			ENDIF
			F= (C*G)+(S*Y)
			X=-(S*G)+(C*Y)
			DO NM=1,M
				Y=A(NM,J)
				Z=A(NM,I)
				A(NM,J)= (Y*C)+(Z*S)
				A(NM,I)=-(Y*S)+(Z*C)
			end do
		end do
		RV1(L)=0.0
		RV1(K)=F
		W(K)=X
	end do
3   CONTINUE
end do
RETURN
END

