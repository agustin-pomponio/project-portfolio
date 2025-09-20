PROGRAM TP_Final 
IMPLICIT NONE

INTEGER op
INTEGER, PARAMETER :: ecuaciones=4
REAL(8), PARAMETER :: pi=3.141592654
REAL(8), PARAMETER :: Re=7000
REAL(8), PARAMETER :: Fpalo=1800 
REAL(8), PARAMETER :: gravedad=9.81
REAL(8) radio, densidad, Cd, area, angulo, alfa, vo_x, vo_y, deltat, masa, paso, vinicial, tmax, roce, tolerancia, ymax
REAL (8) hmin, hmax, hmaxf, hminf
INTEGER (8) iteres, iterem, iterrk4, iterrkf 
REAL(8), DIMENSION(0:ecuaciones) :: v, error

PRINT*, '============================================================================'
PRINT*, '           Trabajo Práctico Final - Análisis Numérico - 2019'
PRINT*, '           Integrantes: Bandiera, Delfina - Malarczuk, Lucila'
PRINT*, '                        Pomponio, Agustin - Sarraute, Martina'
PRINT*, '============================================================================'
PRINT*,

op=6   !Le pongo 6 para que entre la primera iteracion

DO WHILE (op/=0)

tolerancia = 0.000001
paso = 0.1

masa = 0.045
deltat = 0.002
angulo = 20                    
alfa = (angulo*2*pi)/360.      !Como el angulo tiene que estar en radianes para senos y cosenos lo paso

!Componentes de ambas velocidades iniciales, salen de la sumatoria de fuerzas
vo_x = (Fpalo*cos(alfa)*deltat)/masa 
vo_y = (Fpalo*deltat*sin(alfa))/masa
vinicial = sqrt(((vo_x)**2)+(vo_y)**2)

radio = 0.0215
densidad = 1.225
!Calculo del area seccional
area = pi*(radio**2)

!Calculo de Coeficiente de Arrastre
Cd = (((24/Re)**0.52)+(0.32**0.52))**(1/0.52)
PRINT*,
PRINT*, 'Coeficiente de arrastre =', Cd
roce = Cd*0.5*densidad*area

!VALORES INICIALES
!tiempo, posiciones, velocidades
v(0) = 0.001 !Distinto de cero porque t=0 es para el impacto
v(1) = 0.0
v(2) = 0.0
v(3) = vo_x
v(4) = vo_y

PRINT*, 
PRINT*, 'Ingrese el metodo a utilizar'
PRINT*, '0) para salir'
PRINT*, '1) Euler Simple'
PRINT*, '2) Euler Mejorado'
PRINT*, '3) Runge Kutta 4'
PRINT*, '4) Runge Kutta Fehlberg'
READ(*,*) op
CALL SYSTEM ('clear')

tmax=100

SELECT CASE (op)

     CASE(1)

PRINT*, 'Metodo de Euler Simple'
OPEN(UNIT=2, FILE='Trayectoria_ES.dat', STATUS='REPLACE')    !Hay un archivo para cada metodo
WRITE(2,'(2F10.4)') v(1), v(2)
iteres=1
hmax=0                                                       !Para buscar h maximo y minimo
hmin=100
ymax=0
DO WHILE ((v(0)<=tmax) .and. (v(2)>=0.000))                  !Mientras que la pelota no toque el suelo (v(2) es posicion en y)
paso = cambiopaso(v,paso,tolerancia,roce,masa,1,hmax,hmin)   !Si quiero hacer cambio de paso, si lo quiero fijo comento esta linea
v = euler_simple(v,paso,roce,masa)                           !Le tengo que pasar roce y masa porque el metodo llama adentro a vprima (que los necesita)
IF (ymax<v(2)) THEN    !Busca la altura maxima
ymax = v(2)
END IF
WRITE(2,'(2F10.4)') v(1), v(2) 
iteres = iteres + 1
END DO

CLOSE(2)
CALL SYSTEM ("gnuplot -persist script_ES.p")                 !Hay un script para cada metodo

     CASE(2)

PRINT*, 'Metodo de Euler Modificado'
OPEN(UNIT=3, FILE='Trayectoria_EM.dat', STATUS='REPLACE')
WRITE(3,'(2F10.4)') v(1), v(2)
iterem=1
hmax=0
hmin=100
ymax=0
DO WHILE ((v(0)<=tmax) .and. (v(2)>=0.000)) 
paso = cambiopaso(v,paso,tolerancia,roce,masa,2,hmax,hmin)
v = euler_modificado(v,paso,roce,masa)
IF (ymax<v(2)) THEN
ymax = v(2)
END IF
WRITE(3,'(2F10.4)') v(1), v(2)
iterem = iterem + 1
END DO

CLOSE(3)
CALL SYSTEM ("gnuplot -persist script_EM.p")

     CASE(3)

PRINT*, 'Metodo de Runge Kutta 4'
OPEN(UNIT=4, FILE='Trayectoria_RK4.dat', STATUS='REPLACE')
WRITE(4,'(2F10.4)') v(1), v(2)
iterrk4=1
hmax=0
hmin=100
ymax=0
DO WHILE ((v(0)<=tmax) .and. (v(2)>=0.000))
paso = cambiopaso(v,paso,tolerancia,roce,masa,3,hmax,hmin)
v = rk4(v,paso,roce,masa)
IF (ymax<v(2)) THEN
ymax = v(2)
END IF
WRITE(4,'(2F10.4)') v(1), v(2)
iterrk4 = iterrk4 + 1
END DO

CLOSE(4)
CALL SYSTEM ("gnuplot -persist script_RK4.p")

     CASE(4)

error=0.0
PRINT*, 'Metodo de Runge Kutta Fehlberg'
OPEN(UNIT=1, FILE='Trayectoria_RKF.dat', STATUS='REPLACE')
WRITE(1,'(2F10.4)') v(1), v(2)
iterrkf=1
hmaxf=0
hminf=100
DO WHILE ((v(0)<=tmax) .and. (v(2)>=0.000)) 
paso = cambiopaso_rkf(v,paso,tolerancia,roce,masa,error,hmaxf,hminf)
v = rkf(v,paso,roce,masa,error)
IF (ymax<v(2)) THEN    
ymax = v(2)
END IF
WRITE(1,'(2F10.4)') v(1), v(2)
iterrkf = iterrkf + 1
END DO

CLOSE(1)
CALL SYSTEM ("gnuplot -persist script_RKF.p")

END SELECT

PRINT*, 'Altura maxima alcanzada:', ymax
PRINT*, 'Distancia del disparo:', v(1)

END DO

!~ CALL SYSTEM ("gnuplot -persist script_todos.p") !Grafica comparativamente metodos
!~ CALL SYSTEM ("gnuplot -persist script_fuerzas.p") !Grafica comparativamente fuerzas
!~ CALL SYSTEM ("gnuplot -persist script_angulos.p") !Grafica comparativamente angulos


!Para comparacion de angulos: un archivo para cada angulo y despues un script que grafique todos juntos
!Para comparacion de fuerzas: un archivo para cada fuerza y despues un script que grafique todos juntos


CONTAINS

FUNCTION vprima (v,roce,masa)
REAL(8), DIMENSION(0:ecuaciones) :: v, vprima
REAL(8) roce,masa
vprima(0) = 1
vprima(1) = v(3)
vprima(2) = v(4)
vprima(3) = -(roce*(v(3)**2))/masa
vprima(4) = -(roce*(v(4)**2))/masa - gravedad
END FUNCTION

FUNCTION norma(v)
INTEGER i
REAL(8) :: norma
REAL(8), DIMENSION(0:ecuaciones) :: v
norma = ABS(v(0))
DO i = 1, ecuaciones
norma = MAX(norma, ABS(v(i)))
END DO
END FUNCTION

FUNCTION euler_simple (v,h,roce,masa)
REAL(8), INTENT(IN) :: h,roce,masa
REAL(8), DIMENSION(0:ecuaciones) :: v, euler_simple
    euler_simple = v + (h*vprima(v,roce,masa))
END FUNCTION

FUNCTION euler_modificado (v,h,roce,masa)
REAL(8), INTENT(IN) :: h,roce,masa
REAL(8), DIMENSION(0:ecuaciones) :: v, vp, euler_modificado
    vp = vprima(v,roce,masa)
    euler_modificado = v + h*(vp+vprima(v+h*vp,roce,masa))/2.0
END FUNCTION

FUNCTION rk4 (v,h,roce,masa)
REAL(8), INTENT(IN) :: h,roce,masa
REAL(8), DIMENSION(0:ecuaciones) :: rk4, v, k1, k2, k3, k4
    k1 = h*vprima(v,roce,masa)
    k2 = h*vprima(v+k1/2.0,roce,masa)
    k3 = h*vprima(v+k2/2.0,roce,masa)
    k4 = h*vprima(v+k3,roce,masa)
    rk4 = v + (k1 + 2.0*k2 + 2.0*k3 + k4) /6.0
END FUNCTION

FUNCTION rkf (v,h,roce,masa,error)
    REAL(8), INTENT(IN) :: h, roce, masa
    REAL(8), DIMENSION(0:ecuaciones) :: v, k1, k2, k3, k4, k5, k6, error, rkf
    k1 = h * vprima(v,roce,masa)
    k2 = h * vprima(v + k1/4.0, roce, masa)
    k3 = h * vprima(v + (3.0*k1 + 9.0*k2)/32.0, roce, masa)
    k4 = h * vprima(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0, roce, masa)
    k5 = h * vprima(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0, roce, masa)
    k6 = h * vprima(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0, roce, masa)
    rkf = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
    error = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
END FUNCTION


FUNCTION cambiopaso (v,h,tolerancia,roce,masa,metodo,hmax,hmin)
    LOGICAL h_encontrado, h_duplicado, h_dividido
    INTEGER metodo
    REAL(8), INTENT(IN) :: tolerancia, roce, masa 
    REAL(8) :: cambiopaso, error_ajuste, h, hmax, hmin
    REAL(8), DIMENSION(0:ecuaciones) :: v, a, b, c
    
    h_encontrado = .FALSE.
    DO WHILE (.NOT. h_encontrado)  !Itera mientras que no encuentre el h que cumpla tolerancia
    
    IF (metodo==1) THEN   !El metodo lo paso como un numero desde el programa (ES, EM y RK4)
    a = euler_simple(v, h, roce, masa) !Metodo con un paso h
    b = euler_simple(v, h/2, roce, masa) !Metodo con un paso h/2
    c = euler_simple(b, h/2, roce, masa) !Metodo con otro paso h/2
    ELSE IF (metodo==2) THEN
    a = euler_modificado(v, h, roce, masa)
    b = euler_modificado(v, h/2, roce, masa)
    c = euler_modificado(b, h/2, roce, masa)
    ELSE IF (metodo==3) THEN
    a = rk4(v, h, roce, masa)
    b = rk4(v, h/2, roce, masa)
    c = rk4(b, h/2, roce, masa)
    END IF
      
    error_ajuste = norma(a-c)    !Calcula el error como diferencia de pasos con h y h/2
    IF ((error_ajuste > tolerancia) .AND. (.NOT. h_duplicado)) THEN
    h = h / 2
    h_dividido = .TRUE.
    ELSE IF ((error_ajuste < (tolerancia / 2)) .AND. (.NOT. h_dividido)) THEN
    h = h * 2
    h_duplicado = .TRUE.
    ELSE
    h_encontrado = .TRUE.   !Cuando lo encuentra transforma la variable booleana en verdadera
    END IF
    END DO

    cambiopaso = h

    IF (h>hmax) THEN   !Encuentra valor de h minimo y maximo del metodo
        hmax=h
    END IF
    IF(h<hmin) THEN
        hmin=h
    END IF

END FUNCTION


FUNCTION cambiopaso_rkf (v, h, tolerancia, roce, masa, error, hmaxf, hminf)
    REAL(8), INTENT(IN) :: tolerancia, roce, masa
    REAL(8)  cambiopaso_rkf, error_ajuste, h, alfa, hmaxf, hminf
    REAL(8), DIMENSION(0:ecuaciones) :: v, vprox, error

    vprox = rkf(v, h, roce, masa, error) !Hace RKF y le devuelve el vector v y el vector de error 
    error_ajuste = norma(error)
    IF (error_ajuste <= tolerancia) THEN
    alfa = 0.2
    ELSE
    alfa = 0.22
    END IF
    h = h * ((tolerancia / error_ajuste) ** alfa)

    cambiopaso_rkf = h
    
    IF (h<hminf) THEN
        hminf=h
    END IF
    IF (h>hmaxf) THEN
        hmaxf=h
    END IF
    
    
END FUNCTION


END PROGRAM 
