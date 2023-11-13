#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
// OPTIONAL: comment this out for console output
// #define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

#define START_TIMEMEASUREMENT(name)                                                                                         \
  struct timeval __FILE__##__func__##name##actualtime;                                                                      \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);                                                                \
  double __FILE__##__func__##name##s_time = (double)__FILE__##__func__##name##actualtime.tv_sec +                           \
                                            ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0)

#define END_TIMEMEASUREMENT(name, res)                                                                                      \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);                                                                \
  res = (double)__FILE__##__func__##name##actualtime.tv_sec +                                                               \
        ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0) - __FILE__##__func__##name##s_time

typedef uint8_t number_type;
typedef uint64_t header_type;
#define NUMBER_TYPE_VTK_NAME "UInt8"
#define HEADER_TYPE_VTK_NAME "UInt64"

void myexit(const char *s, ...) {
  va_list args;
  va_start(args, s);
  vprintf(s, args);
  printf("\n");
  va_end(args);
  abort();
}

int testLittleEndian() {
  int32_t test = 1;
  char *testdata = (char *)&test;
  if (testdata[0] == 1) {
    return 1;
  }
  return 0;
}

const char *vtk_header_template =
    "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"%s\" header_type=\"" HEADER_TYPE_VTK_NAME "\">\n"
    "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 0\">\n"
    "    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
    "      <PointData>\n"
    "      </PointData>\n"
    "      <CellData Scalars=\"GameOfLife\">\n"
    "        <DataArray type=\"" NUMBER_TYPE_VTK_NAME
    "\" Name=\"GameOfLife\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n"
    "      </CellData>\n"
    "    </Piece>\n"
    "  </ImageData>\n"
    "  <AppendedData encoding=\"raw\">\n"
    "   _";

char vtk_header[10000];
void create_vtk_header(char *header, int width, int height) {
  snprintf(header, 10000, vtk_header_template, testLittleEndian() ? "LittleEndian" : "BigEndian", width, height, width,
           height);
}

char *vtk_tail = "\n  </AppendedData>\n"
                 "</VTKFile>\n";

void write_vtk_data(FILE *f, char *data, int length) {
  if (fwrite(data, 1, length, f) != length) {
    myexit("Could not write vtk-Data");
  }
}

void write_field(number_type *currentfield, int width, int height, int timestep) {
#ifdef CONSOLE_OUTPUT
  printf("\033[H");
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++)
      printf(ALIVE == currentfield[calcIndex(width, x, y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
  printf("\ntimestep=%d", timestep);
  usleep(80000);
#else
  if (timestep == 0) {
    mkdir("./gol/", 0777);
    create_vtk_header(vtk_header, width, height);
  }
  printf("writing timestep %d\n", timestep);
  FILE *fp; // The current file handle.
  char filename[1024];
  snprintf(filename, 1024, "./gol/gol-%05d.vti", timestep);
  fp = fopen(filename, "w");
  // write the header of the vti file
  write_vtk_data(fp, vtk_header, strlen(vtk_header));
  // write the length in bytes before writing the data
  header_type length = width * height * sizeof(number_type);
  write_vtk_data(fp, (char *)&length, sizeof(header_type));
  write_vtk_data(fp, (char *)currentfield, width * height * sizeof(number_type));
  // write the tail of the vti file
  write_vtk_data(fp, vtk_tail, strlen(vtk_tail));
  fclose(fp);
  printf("finished writing timestep %d\n", timestep);
#endif
}


void evolve(number_type *currentfield, number_type *newfield, int width, int height) {
  /**
   * navigate to current cell, remember its index, and start counter for living neighbours
  */
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      int currentCellIndex = calcIndex(width, x, y);
      int noOfLivingNeighbours = 0;

      // count living neighbours around current cell
      /** neighbour cell layout:
       * (-1,-1)   (-1,0)    (-1,1)
       * (0,-1)     (0,0)     (0,1)   <-- (0,0) is active cell, needs to be skipped
       * (1,-1)     (1,0)     (1,1)
      */
      for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
          if (dx == 0 && dy == 0) continue; // skip the active cell
          if (currentfield[calcIndex(width, x + dx, y + dy)] == ALIVE) {
            noOfLivingNeighbours++;
          }
        }
      }

      /** if necessary, update state of current field according to rules:
       * ALIVE && < 2 neighbours -> DEAD
       * ALIVE && 2 or 3 neighbours -> stays ALIVE
       * ALIVE && > 3 neighbours -> DEAD
       * DEAD && == 3 neighbours -> ALIVE
      */
      newfield[currentCellIndex] = (currentfield[currentCellIndex] == ALIVE && (noOfLivingNeighbours == 2 || noOfLivingNeighbours == 3)) || 
                        (currentfield[currentCellIndex] == DEAD && noOfLivingNeighbours == 3) ? ALIVE : DEAD;
    }
  }
}

void filling_random(number_type *currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner(number_type *currentfield, int width, int height) {
  int offset_x = width / 3;
  int offset_y = height / 2;
  currentfield[calcIndex(width, offset_x + 0, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 1, offset_y + 2)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 0)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 2)] = ALIVE;
}

void apply_periodic_boundaries(number_type *field, int width, int height) {
  /**
   * the new field does not yet have the boundary logic implemented
   * copies edges to all four sides of the field to implement periodic boundaries:
  */
    for (int y = 1; y < height - 1; y++) {
        // copy left to right and right to left
        field[calcIndex(width, 0, y)] = field[calcIndex(width, width - 2, y)];
        field[calcIndex(width, width - 1, y)] = field[calcIndex(width, 1, y)];
    }

      for (int x = 0; x < width; x++) {
        // copy top to bottom and bottom to top
        field[calcIndex(width, x, 0)] = field[calcIndex(width, x, height - 2)];
        field[calcIndex(width, x, height - 1)] = field[calcIndex(width, x, 1)];
    }
}

void game(int width, int height, int num_timesteps) {
  number_type *currentfield = calloc(width * height, sizeof(number_type));
  number_type *newfield = calloc(width * height, sizeof(number_type));

  // filling_random (currentfield, width, height);
  filling_runner(currentfield, width, height);

  int time = 0;
  write_field(currentfield, width, height, time);
  apply_periodic_boundaries(currentfield, width, height);

  for (time = 1; time <= num_timesteps; time++) {
    evolve(currentfield, newfield, width, height);
    write_field(newfield, width, height, time);
    apply_periodic_boundaries(newfield, width, height);

    // swap fields by swapping the pointers
    number_type *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v) {
  int width = 0, height = 0, num_timesteps;
  if (c == 4) {
    width = atoi(v[1]) + 2;     ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]) + 2;    ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps
    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }

    double elapsed_time;
    START_TIMEMEASUREMENT(measure_game_time);

    game(width, height, num_timesteps);

    END_TIMEMEASUREMENT(measure_game_time, elapsed_time);
    printf("time elapsed: %lf sec\n", elapsed_time);
  } else {
    myexit("Too few arguments, example: ./gameoflife <x size> <y size> <number of timesteps>");
  }
}
