#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include <ctype.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define BLOCK_SIZE 0.5f
#define CHUNK_SIZE 16
#define CHUNK_HEIGHT 128 
#define WORLD_RADIUS 4
#define MAX_CHUNKS 1024
#define MAX_BLOCKS_PER_CHUNK (CHUNK_SIZE*CHUNK_SIZE*CHUNK_HEIGHT)
#define PLAYER_HEIGHT (BLOCK_SIZE * 3.5f) 
#define MAX_STEP_HEIGHT (BLOCK_SIZE * 0.5f)
#define GRAVITY -20.0f
#define BLOCKS_PER_SECOND 7.0f 
#define MOVE_SPEED (BLOCKS_PER_SECOND * BLOCK_SIZE)  
#define FLY_SPEED (BLOCKS_PER_SECOND * BLOCK_SIZE * 2.0f) 
#define JUMP_HEIGHT_BLOCKS 2.0f
#define JUMP_SPEED 7.0f 
#define DAY_LENGTH 60.0f 
typedef struct { float x,y,z,r,g,b; } Block;
typedef int BlockMap[CHUNK_SIZE][CHUNK_HEIGHT][CHUNK_SIZE]; 
typedef struct { 
    int cx, cz; 
    int count; 
    Block* blocks; 
    BlockMap map;  
    GLuint dlist; 
} Chunk;
static Chunk chunks[MAX_CHUNKS];
static int chunk_count=0;
static float player_x=0,player_y=0,player_z=0;
static float vel_x=0,vel_y=0,vel_z=0;
static float yaw=180,pitch=-6;
static int key_w=0,key_s=0,key_a=0,key_d=0,key_space=0,key_c=0;
static int is_flying=0,on_ground=0;
static int win_w=1280,win_h=720;
static float world_time = 0.0f; 
void generate_tree(Chunk* c, int lx, int ly, int lz);
void generate_house(Chunk* c, int lx, int ly, int lz);
void generate_flower(Chunk* c, int lx, int ly, int lz);
unsigned hash(unsigned x){ x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;x^=x>>16;return x;}
unsigned hash_coords(int x,int z){ return hash((unsigned)x*374761393u ^ (unsigned)z*668265263u ^ 0x9e3779b1u);}
float hash01(unsigned v){ return (v%10000)/10000.0f;}
float smoothstep(float t){ return t*t*(3-2*t);}
float noise(float x,float z){
    int ix=(int)floor(x), iz=(int)floor(z);
    float tx=x-ix,tz=z-iz;
    unsigned h1=hash_coords(ix,iz),h2=hash_coords(ix+1,iz),h3=hash_coords(ix,iz+1),h4=hash_coords(ix+1,iz+1);
    float f1=hash01(h1),f2=hash01(h2),f3=hash01(h3),f4=hash01(h4);
    float sx=smoothstep(tx);
    float a=f1+sx*(f2-f1), b=f3+sx*(f4-f3), sy=smoothstep(tz);
    return a + sy*(b-a);
}
float get_ground_height(float x,float z){
    float h=0;
    h+=noise(x*0.005f,z*0.005f)*30; 
    h+=noise(x*0.02f,z*0.02f)*10;  
    h+=noise(x*0.1f,z*0.1f)*3;     
    if(h<1) h=1; if(h>80) h=80;
    return h;
}
Chunk* find_chunk(int cx,int cz){
    for(int i=0;i<chunk_count;i++) if(chunks[i].cx==cx && chunks[i].cz==cz) return &chunks[i];
    return NULL;
}
Chunk* make_chunk(int cx,int cz){
    if(chunk_count>=MAX_CHUNKS){ fprintf(stderr,"Too many chunks\n"); exit(1); }
    Chunk* c=&chunks[chunk_count++];
    c->cx=cx; c->cz=cz; c->count=0;
    c->blocks=malloc(sizeof(Block)*MAX_BLOCKS_PER_CHUNK);
    for(int x=0; x<CHUNK_SIZE; x++) for(int y=0; y<CHUNK_HEIGHT; y++) for(int z=0; z<CHUNK_SIZE; z++) c->map[x][y][z] = 0;
    c->dlist=0;
    return c;
}
void add_block_to_map(Chunk* c, int lx, int ly, int lz, float r, float g, float b){
    if(lx>=0 && lx<CHUNK_SIZE && ly>=0 && ly<CHUNK_HEIGHT && lz>=0 && lz<CHUNK_SIZE && c->count<MAX_BLOCKS_PER_CHUNK){
        float wx_corner = (c->cx*CHUNK_SIZE + lx) * BLOCK_SIZE;
        float wz_corner = (c->cz*CHUNK_SIZE + lz) * BLOCK_SIZE;
        float wy_corner = ly * BLOCK_SIZE;
        float wx = wx_corner + BLOCK_SIZE / 2.0f;
        float wy = wy_corner + BLOCK_SIZE / 2.0f;
        float wz = wz_corner + BLOCK_SIZE / 2.0f;
        c->blocks[c->count++] = (Block){wx,wy,wz,r,g,b};
        c->map[lx][ly][lz] = 1; 
    }
}
#define IS_SOLID(nx, ny, nz) (nx>=0 && nx<CHUNK_SIZE && ny>=0 && ny<CHUNK_HEIGHT && nz>=0 && nz<CHUNK_SIZE && c->map[nx][ny][nz])
void draw_face(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, float x4, float y4, float z4) {
    glBegin(GL_QUADS);
        glVertex3f(x1, y1, z1);
        glVertex3f(x2, y2, z2);
        glVertex3f(x3, y3, z3);
        glVertex3f(x4, y4, z4);
    glEnd();
}
void draw_block_optimized(const Block* b, int lx, int ly, int lz, const Chunk* c){
    glColor3f(b->r, b->g, b->b);
    float X = b->x, Y = b->y, Z = b->z;
    float S = BLOCK_SIZE;
    if (!IS_SOLID(lx, ly + 1, lz)) {
        glNormal3f(0, 1, 0);
        draw_face(X-S/2, Y+S/2, Z-S/2, X+S/2, Y+S/2, Z-S/2, X+S/2, Y+S/2, Z+S/2, X-S/2, Y+S/2, Z+S/2);
    }
    if (!IS_SOLID(lx, ly - 1, lz)) {
        glNormal3f(0, -1, 0);
        draw_face(X-S/2, Y-S/2, Z-S/2, X-S/2, Y-S/2, Z+S/2, X+S/2, Y-S/2, Z+S/2, X+S/2, Y-S/2, Z-S/2);
    }
    if (!IS_SOLID(lx, ly, lz + 1)) {
        glNormal3f(0, 0, 1);
        draw_face(X-S/2, Y-S/2, Z+S/2, X+S/2, Y-S/2, Z+S/2, X+S/2, Y+S/2, Z+S/2, X-S/2, Y+S/2, Z+S/2);
    }
    if (!IS_SOLID(lx, ly, lz - 1)) {
        glNormal3f(0, 0, -1);
        draw_face(X-S/2, Y-S/2, Z-S/2, X-S/2, Y+S/2, Z-S/2, X+S/2, Y+S/2, Z-S/2, X+S/2, Y-S/2, Z-S/2);
    }
    if (!IS_SOLID(lx + 1, ly, lz)) {
        glNormal3f(1, 0, 0);
        draw_face(X+S/2, Y-S/2, Z-S/2, X+S/2, Y+S/2, Z-S/2, X+S/2, Y+S/2, Z+S/2, X+S/2, Y-S/2, Z+S/2);
    }
    if (!IS_SOLID(lx - 1, ly, lz)) {
        glNormal3f(-1, 0, 0);
        draw_face(X-S/2, Y-S/2, Z-S/2, X-S/2, Y-S/2, Z+S/2, X-S/2, Y+S/2, Z+S/2, X-S/2, Y+S/2, Z-S/2);
    }

    #undef IS_SOLID
}
void generate_tree(Chunk* c, int lx, int ly, int lz){
    int h=3+(hash_coords((c->cx*CHUNK_SIZE+lx),(c->cz*CHUNK_SIZE+lz))%3);
    for(int i=0;i<h;i++) add_block_to_map(c,lx,ly+i,lz,0.55f,0.27f,0.07f);
    for(int dx=-2;dx<=2;dx++) for(int dz=-2;dz<=2;dz++) for(int dy=0;dy<=2;dy++){
        if(fabs(dx)+fabs(dz)+fabs(dy)<4) add_block_to_map(c,lx+dx,ly+h+dy,lz+dz,0.12f,0.6f,0.12f);
    }
}
void generate_flower(Chunk* c, int lx, int ly, int lz){
    int id=hash_coords((c->cx*CHUNK_SIZE+lx),(c->cz*CHUNK_SIZE+lz))%3;
    if(id==0) add_block_to_map(c,lx,ly,lz,1,0.2f,0.2f);
    else if(id==1) add_block_to_map(c,lx,ly,lz,1,1,0.2f);
    else add_block_to_map(c,lx,ly,lz,0.4f,0.4f,1);
}
void generate_house(Chunk* c, int lx, int ly, int lz){
    int w=3+(hash_coords((c->cx*CHUNK_SIZE+lx+1),(c->cz*CHUNK_SIZE+lz))%3);
    int l=3+(hash_coords((c->cx*CHUNK_SIZE+lx+2),(c->cz*CHUNK_SIZE+lz))%3);
    int h=2+(hash_coords((c->cx*CHUNK_SIZE+lx+3),(c->cz*CHUNK_SIZE+lz))%2);
    for(int px=0;px<w;px++) for(int pz=0;pz<l;pz++) for(int py=0;py<h;py++){
        if(px==0||px==w-1||pz==0||pz==l-1||py==0||py==h-1) add_block_to_map(c,lx+px,ly+py,lz+pz,0.85f,0.75f,0.6f);
    }
    for(int px=-1;px<w+1;px++) for(int pz=-1;pz<l+1;pz++) add_block_to_map(c,lx+px,ly+h,lz+pz,0.6f,0.2f,0.2f);
}
void generate_chunk(int cx,int cz){
    if(find_chunk(cx,cz)) return;
    Chunk* c=make_chunk(cx,cz);
    for(int x=0;x<CHUNK_SIZE;x++) for(int z=0;z<CHUNK_SIZE;z++){
        float wx_sample=(cx*CHUNK_SIZE+x)*BLOCK_SIZE;
        float wz_sample=(cz*CHUNK_SIZE+z)*BLOCK_SIZE;
        float h=get_ground_height(wx_sample,wz_sample);
        
        int ih=(int)floor(h/BLOCK_SIZE); 
        
        for(int y=0;y<=ih;y++){
            float r,g,b;
            if(y==ih) { r=0.2f; g=0.8f; b=0.2f; }
            else if(y>ih-4) { r=0.6f; g=0.4f; b=0.3f; }
            else { r=0.5f; g=0.5f; b=0.5f; }
            
            add_block_to_map(c,x,y,z,r,g,b);
        }
        unsigned hv=hash_coords((int)wx_sample,(int)wz_sample);
        int sy = ih + 1;
        if(sy < CHUNK_HEIGHT){
            if(hv%40==0) generate_tree(c,x,sy,z);
            else if(hv%60==0) generate_house(c,x,sy,z);
            else if(hv%20==0) generate_flower(c,x,sy,z);
        }
    }
    c->dlist=glGenLists(1);
    glNewList(c->dlist,GL_COMPILE);
    for(int i=0;i<c->count;i++){
        int lx = (int)roundf(c->blocks[i].x/BLOCK_SIZE - 0.5f) - c->cx*CHUNK_SIZE;
        int lz = (int)roundf(c->blocks[i].z/BLOCK_SIZE - 0.5f) - c->cz*CHUNK_SIZE;
        int ly = (int)roundf(c->blocks[i].y/BLOCK_SIZE - 0.5f);
        
        draw_block_optimized(&c->blocks[i], lx, ly, lz, c);
    }
    glEndList();
    
    free(c->blocks);
    c->blocks = NULL;
}
void ensure_chunks(float px,float pz){
    int pcx=(int)floor(px/(CHUNK_SIZE*BLOCK_SIZE)), pcz=(int)floor(pz/(CHUNK_SIZE*BLOCK_SIZE));
    for(int dx=-WORLD_RADIUS;dx<=WORLD_RADIUS;dx++) for(int dz=-WORLD_RADIUS;dz<=WORLD_RADIUS;dz++)
        if(!find_chunk(pcx+dx,pcz+dz)) generate_chunk(pcx+dx,pcx+dz); // FIX: Should use dz here, not pcx+dz
}
void draw_world(){ for(int i=0;i<chunk_count;i++) glCallList(chunks[i].dlist); }
void ensure_chunks_fixed(float px,float pz){
    int pcx=(int)floor(px/(CHUNK_SIZE*BLOCK_SIZE)), pcz=(int)floor(pz/(CHUNK_SIZE*BLOCK_SIZE));
    for(int dx=-WORLD_RADIUS;dx<=WORLD_RADIUS;dx++) for(int dz=-WORLD_RADIUS;dz<=WORLD_RADIUS;dz++)
        if(!find_chunk(pcx+dx,pcz+dz)) generate_chunk(pcx+dx,pcz+dz);
}


void update_lighting(float dt) {
    world_time += dt;
    if (world_time > DAY_LENGTH) world_time -= DAY_LENGTH;
    float time_factor = sinf(world_time / DAY_LENGTH * 2.0f * M_PI);
    time_factor = (time_factor * 0.5f) + 0.5f; 
    float min_light = 0.2f;
    float max_light = 1.0f;
    float light_intensity = min_light + (max_light - min_light) * time_factor;
    float sky_day[] = {0.5f, 0.7f, 1.0f};
    float sky_night[] = {0.05f, 0.05f, 0.15f};
    float sky_r = sky_night[0] + (sky_day[0] - sky_night[0]) * time_factor;
    float sky_g = sky_night[1] + (sky_day[1] - sky_night[1]) * time_factor;
    float sky_b = sky_night[2] + (sky_day[2] - sky_night[2]) * time_factor;
    glClearColor(sky_r, sky_g, sky_b, 1.0f);
    GLfloat ambient[] = { light_intensity * 0.7f, light_intensity * 0.7f, light_intensity * 0.7f, 1.0f };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    GLfloat diffuse[] = { light_intensity, light_intensity, light_intensity, 1.0f };
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    float angle = world_time / DAY_LENGTH * 2.0f * M_PI;
    GLfloat light_pos[] = { sinf(angle) * 50.0f, cosf(angle) * 50.0f, 0.0f, 0.0f }; 
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
}


void update_player(float dt){
    update_lighting(dt);
    float fx=0,fz=0;
    float speed = is_flying ? FLY_SPEED : MOVE_SPEED;
    if(key_w) fx+=sinf(yaw*M_PI/180.0f), fz+=cosf(yaw*M_PI/180.0f);
    if(key_s) fx-=sinf(yaw*M_PI/180.0f), fz-=cosf(yaw*M_PI/180.0f);
    if(key_a) fx-=cosf(yaw*M_PI/180.0f), fz+=sinf(yaw*M_PI/180.0f);
    if(key_d) fx+=cosf(yaw*M_PI/180.0f), fz-=sinf(yaw*M_PI/180.0f);
    float next_x = player_x + fx*speed*dt;
    float next_z = player_z + fz*speed*dt;
    float current_ground_y = get_ground_height(player_x, player_z) + PLAYER_HEIGHT;
    float new_ground_x = get_ground_height(next_x, player_z) + PLAYER_HEIGHT;
    float new_ground_z = get_ground_height(player_x, next_z) + PLAYER_HEIGHT;
    if (new_ground_x <= current_ground_y + MAX_STEP_HEIGHT) {
        player_x = next_x;
        if (!is_flying && new_ground_x > current_ground_y && on_ground) {
             player_y = new_ground_x;
        }
    }
    if (new_ground_z <= current_ground_y + MAX_STEP_HEIGHT) {
        player_z = next_z;
        if (!is_flying && new_ground_z > current_ground_y && on_ground) {
             player_y = new_ground_z;
        }
    }
    if(is_flying){
        if(key_space) player_y+=FLY_SPEED*dt;
        if(key_c) player_y-=FLY_SPEED*dt;
        vel_y=0; on_ground=0;
    } else {
        vel_y+=GRAVITY*dt;
        if(on_ground && key_space){ vel_y=JUMP_SPEED; on_ground=0; }
        player_y+=vel_y*dt;
        
        float ground=get_ground_height(player_x,player_z)+PLAYER_HEIGHT;
        if(player_y<ground){ 
            player_y=ground; 
            vel_y=0; 
            on_ground=1; 
        } else {
            on_ground = 0;
        }
    }
    ensure_chunks_fixed(player_x,player_z);
}
void display(){
    static int last_time=0;
    int time_ms=glutGet(GLUT_ELAPSED_TIME);
    float dt=(time_ms-last_time)/1000.0f;
    last_time=time_ms;
    update_player(dt);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    float yaw_rad=yaw*M_PI/180.0f, pitch_rad=pitch*M_PI/180.0f;
    float dx=cosf(pitch_rad)*sinf(yaw_rad), dy=sinf(pitch_rad), dz=cosf(pitch_rad)*cosf(yaw_rad);
    gluLookAt(player_x,player_y,player_z, player_x+dx,player_y+dy,player_z+dz,0,1,0);
    draw_world();
    glutSwapBuffers();
    glutPostRedisplay();
}
void reshape(int w,int h){ 
    win_w=w; win_h=h;
    glViewport(0,0,w,h); 
    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity(); 
    gluPerspective(65,(float)w/h,0.1,1000); 
    glMatrixMode(GL_MODELVIEW); 
    glutWarpPointer(win_w / 2, win_h / 2);
}
void key_down(unsigned char k,int x,int y){ 
    k=tolower(k); 
    if(k=='w') key_w=1; 
    if(k=='s') key_s=1; 
    if(k=='a') key_a=1; 
    if(k=='d') key_d=1; 
    if(k==' ') key_space=1; 
    if(k=='c') key_c=1; 
    if(k=='f') is_flying=!is_flying;
    if(k==27) exit(0);
}
void key_up(unsigned char k,int x,int y){ 
    k=tolower(k); 
    if(k=='w') key_w=0; 
    if(k=='s') key_s=0; 
    if(k=='a') key_a=0; 
    if(k=='d') key_d=0; 
    if(k==' ') key_space=0; 
    if(k=='c') key_c=0;
}
void mouse_motion(int x,int y){ 
    int dx = x - win_w / 2;
    int dy = y - win_h / 2;
    if (dx != 0 || dy != 0) {
        yaw += dx * 0.2f;
        pitch -= dy * 0.2f;
        if(pitch > 89) pitch = 89;
        if(pitch < -89) pitch = -89;
        glutWarpPointer(win_w / 2, win_h / 2);
    }
}
int main(int argc,char**argv){
    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    glutInitWindowSize(win_w,win_h);
    glutCreateWindow("minceraft");
    glEnable(GL_DEPTH_TEST); 
    glEnable(GL_LIGHTING); 
    glEnable(GL_LIGHT0); 
    glEnable(GL_COLOR_MATERIAL); 
    glShadeModel(GL_SMOOTH); 
    update_lighting(0.0f); 
    glutSetCursor(GLUT_CURSOR_NONE);
    glutWarpPointer(win_w/2, win_h/2);
    ensure_chunks_fixed(player_x,player_z);
    float initial_height = get_ground_height(player_x,player_z);
    int initial_block_y = (int)floor(initial_height / BLOCK_SIZE);
    player_y = initial_block_y * BLOCK_SIZE + PLAYER_HEIGHT + 0.1f;
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key_down);
    glutKeyboardUpFunc(key_up);
    glutPassiveMotionFunc(mouse_motion);
    glutMainLoop();
    for(int i=0; i<chunk_count; i++) {
        if(chunks[i].dlist) glDeleteLists(chunks[i].dlist, 1);
    }
    return 0;
}