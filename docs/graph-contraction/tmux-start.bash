#!/bin/sh
SESSION="graph-contraction"
PDFREADER="xdg-open"

tmux -2 new-session -d -s $SESSION

tmux new-window -t $SESSION:1 -k -n main
tmux send-keys -t $SESSION:1 'vim contraction.tex' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe body.tex' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe fig1.tex' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe fig2.tex' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe fig3.tex' C-m
tmux send-keys -t $SESSION:1 ':' 'tabe fig4.tex' C-m
tmux send-keys -t $SESSION:1 'gt' C-m
tmux send-keys -t $SESSION:1 'gt' C-m
tmux split-window -h -p 20
tmux send-keys -t $SESSION:1 $PDFREADER ' contraction.pdf &' C-m
tmux select-pane -t 0

tmux new-window -t $SESSION:2 -n make
tmux select-window -t $SESSION:2
tmux send-keys -t $SESSION:2 'vim tmux-start.bash' C-m
tmux send-keys -t $SESSION:2 ':' 'tabe makefile' C-m
tmux split-window -h
tmux select-pane -t 0

tmux select-window -t $SESSION:1

tmux attach -t $SESSION
