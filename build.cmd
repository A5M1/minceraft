@echo off
title building
gcc miner.c -o dist\demo.exe -lfreeglut -lopengl32 -lglu32  -Os -mwindows
echo built!
title done
pause