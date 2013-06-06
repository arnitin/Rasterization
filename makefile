all :
	gcc glm.c -c -o glm.obj -lGL -lGLU -lglut -nostartfiles -g -lm
	gcc gltb.c -c -o gltb.obj -lGL -lGLU -lglut -nostartfiles -g -lm
	gcc smooth.c -o smooth glm.obj gltb.obj -lGL -lGLU -lglut
clean :
	 rm -rf *o *out smooth