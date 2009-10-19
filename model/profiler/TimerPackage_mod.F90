module TimerPackage_mod
   use Timer_mod, only: Timer_type
   use Timer_mod, only: start, stop, reset
   use Timer_mod, only: summary
   use Timer_mod, only: TIMER_SUMMARY_LENGTH
#ifdef USE_MPI
   use Timer_mod, only: parallelSummary
#endif

   use TimerList_mod, only: addTimer, getTimer
   use TimerList_mod, only: start, stop
   use TimerList_mod, only: initialize, finalize, reset
   use TimerList_mod, only: printSummary
#ifdef USE_MPI
   use TimerList_mod, only: printParallelSummary
#endif

end module TimerPackage_mod
