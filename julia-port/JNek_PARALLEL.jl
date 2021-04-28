#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_PARALLEL

      using MPI

      export nid0,
             rank, 
             comm


        const nid0 = 0

        function __init__()

          global comm, rank

          MPI.Init()

          comm = MPI.COMM_WORLD
          rank = MPI.Comm_rank(comm)

          if MPI.Comm_rank(comm) == nid0
            println("Started MPI with $(MPI.Comm_size(comm)) ranks\n")
          end
          
          MPI.Barrier(comm)

        end     # __init__

      end   # Module JNek_PARALLEL











