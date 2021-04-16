#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_MPI

      using MPI

      export nid0, comm, MPI

      const nid0 = 0

      function __init__()

        global comm

        MPI.Init()

        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Started MPI with $(MPI.Comm_size(comm)) ranks\n")
        end
        
        MPI.Barrier(comm)

      end     # jnek_mpi_init

      end   # Module JNek_MPI











