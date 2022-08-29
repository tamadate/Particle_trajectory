NAME        = ../../testTrajectory/trajectory.out
SRCS        = ./src/*cpp
OBJS        = $(SRCS:.cpp=.o)
CXX         = g++
CXXFLAGS    = -std=c++11

.PHONY: all
all: $(NAME)

$(NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(NAME) $(OBJS)

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: fclean
fclean: clean
	$(RM) $(NAME)

.PHONY: re
re: fclean all
