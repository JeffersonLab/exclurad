OBJ=  exclurad.o fint.o

#CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib  -lmathlib  -lpacklib 
 
exclurad : $(OBJ) $(OBJC)
	mkdir -p build
	gfortran -o build/exclurad $(addprefix build/, $(OBJ)) build/$(OBJC) 
$(OBJ) : %.o: %.F
	mkdir -p build
	gfortran -fno-automatic -ffixed-line-length-none -fno-second-underscore -funroll-loops -fomit-frame-pointer  -c $< -o build/$@  
$(OBJC) : %.o: %.c
	mkdir -p build
	gcc -c $< -o build/$@ 
clean:
	rm -f build/exclurad $(addprefix build/, $(OBJ)) $(addprefix build/, $(OBJC)) 



##g77  -O2 -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore \
##        -DLinux \
##        -I. -I./ -I/group/clas/builds/release-4-14/packages/include -I/group/clas/builds/release-4-14/packages/inc_derived -I/apps/tcl/include -I/usr/X11R6/include -c \
##        aac.F  -o /home/avakian/w6/tmp/obj/LinuxRHEL3/generator/aac.o

