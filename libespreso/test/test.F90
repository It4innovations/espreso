PROGRAM test

  USE FETI4I

  TYPE(C_PTR) :: instance
  TYPE(C_PTR) :: matrix

  INTEGER(KIND=FETI4IInt) :: n = 0
  INTEGER(KIND=FETI4IInt) :: eSize = 0
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: indices(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: values(:)

  INTEGER(KIND=FETI4IInt) :: rhs_size = 0
  INTEGER(KIND=FETI4IInt) :: dirichlet_size = 0
  INTEGER(KIND=FETI4IInt) :: neighbours_size = 0

  REAL(KIND=FETI4IReal), ALLOCATABLE :: rhs(:)
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: l2g(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: solution(:)
  INTEGER(KIND=FETI4IInt), ALLOCATABLE :: dirichlet_indices(:)
  REAL(KIND=FETI4IReal), ALLOCATABLE :: dirichlet_values(:)
  INTEGER(KIND=FETI4IMPIInt), ALLOCATABLE :: neighbours(:)


  CALL MPI_INIT()

  CALL TEST4IGetElementsInfo(n, eSize)

  CALL FETI4ICreateStiffnessMatrix(matrix, 0)
  ALLOCATE(values(n*n), indices(n))
  DO i = 1, n
    CALL TEST4IGetElement(i - 1, indices, values)
    CALL FETI4IAddElement(matrix, eSize, indices, values)
  END DO
  write (*,*) "K created"

  CALL TEST4IGetInstanceInfo(rhs_size, dirichlet_size, neighbours_size)

  ALLOCATE(rhs(rhs_size), l2g(rhs_size), solution(rhs_size))
  ALLOCATE(dirichlet_indices(dirichlet_size), dirichlet_values(dirichlet_size))
  ALLOCATE(neighbours(neighbours_size))

  CALL TEST4IGetInstance(rhs, l2g, dirichlet_indices, dirichlet_values, neighbours)

  CALL FETI4ICreateInstance(instance, matrix, rhs_size, rhs, l2g, &
    neighbours_size, neighbours, dirichlet_size, dirichlet_indices, dirichlet_values)
  write (*,*) "Initialized"

  CALL FETI4ISolve(instance, rhs_size, solution)
  write (*,*) "SOLVER ;-D"

  CALL MPI_FINALIZE()

END PROGRAM test 
