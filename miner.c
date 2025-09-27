#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include <ctype.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BLOCK_SIZE 1.0f
#define CHUNK_SIZE 16
#define CHUNK_HEIGHT 128
#define WORLD_RADIUS 4
#define MAX_CHUNKS 1024
#define MAX_BLOCKS_PER_CHUNK (CHUNK_SIZE*CHUNK_SIZE*CHUNK_HEIGHT)
#define PLAYER_HEIGHT (BLOCK_SIZE * 1.8f)
#define PLAYER_RADIUS 0.3f
#define MAX_STEP_HEIGHT (BLOCK_SIZE * 0.5f)
#define GRAVITY -20.0f
#define BLOCKS_PER_SECOND 7.0f
#define MOVE_SPEED (BLOCKS_PER_SECOND * BLOCK_SIZE)
#define FLY_SPEED (BLOCKS_PER_SECOND * BLOCK_SIZE * 2.0f)
#define JUMP_SPEED 7.0f
#define DAY_LENGTH 180.0f

typedef enum{MAT_AIR=0,MAT_GRASS,MAT_DIRT,MAT_STONE,MAT_WOOD,MAT_LEAVES,MAT_FLOWER_RED,MAT_FLOWER_YELLOW,MAT_FLOWER_BLUE,MAT_BRICK,MAT_ROOF,MAT_ANIMAL}Material;

typedef struct{float r,g,b;int solid;}MaterialInfo;

static MaterialInfo material_table[]={
    {0,0,0,0},
    {0.2f,0.8f,0.2f,1},
    {0.6f,0.4f,0.3f,1},
    {0.5f,0.5f,0.5f,1},
    {0.7f,0.5f,0.3f,1},
    {0.12f,0.6f,0.12f,1},
    {1.0f,0.2f,0.2f,0},
    {1.0f,1.0f,0.2f,0},
    {0.4f,0.4f,1.0f,0},
    {0.6f,0.3f,0.3f,1},
    {0.6f,0.2f,0.2f,1},
    {1.0f,0.5f,0.0f,0}
};

typedef struct{float x,y,z;Material mat;}Block;
typedef int BlockMap[CHUNK_SIZE][CHUNK_HEIGHT][CHUNK_SIZE];
typedef struct{int cx,cz;int count;Block*blocks;BlockMap map;GLuint dlist;}Chunk;

static Chunk chunks[MAX_CHUNKS];
static int chunk_count=0;
static float player_x=0,player_y=0,player_z=0;
static float vel_x=0,vel_y=0,vel_z=0;
static float yaw=180,pitch=-6;
static int key_w=0,key_s=0,key_a=0,key_d=0,key_space=0,key_c=0;
static int is_flying=0,on_ground=0;
static int win_w=1280,win_h=720;
static float world_time=0.0f;
static int paused=0;

float q_rsqrt(float number){
    long i;
    float x2,y;
    const float threehalfs=1.5F;
    y=number*0.5F;
    i=*(long*)&number;
    i=0x5f3759df-(i>>1);
    number=*(float*)&i;
    number=number*(threehalfs-(x2=y*number*number));
    return number;
}

unsigned hash(unsigned x){
    x^=x>>16;
    x*=0x7feb352du;
    x^=x>>15;
    x*=0x846ca68bu;
    x^=x>>16;
    return x;
}

unsigned hash_coords(int x,int z){
    return hash((unsigned)x*374761393u ^ (unsigned)z*668265263u ^ 0x9e3779b1u);
}

float hash01(unsigned v){
    return (v%10000)/10000.0f;
}

float smoothstep(float t){
    return t*t*(3-2*t);
}

float noise(float x,float z){
    int ix=(int)floor(x),iz=(int)floor(z);
    float tx=x-ix,tz=z-iz;
    unsigned h1=hash_coords(ix,iz),h2=hash_coords(ix+1,iz),h3=hash_coords(ix,iz+1),h4=hash_coords(ix+1,iz+1);
    float f1=hash01(h1),f2=hash01(h2),f3=hash01(h3),f4=hash01(h4);
    float sx=smoothstep(tx);
    float a=f1+sx*(f2-f1),b=f3+sx*(f4-f3),sy=smoothstep(tz);
    return a+sy*(b-a);
}

float get_ground_height(float x,float z){
    float h=0;
    h+=noise(x*0.005f,z*0.005f)*30;
    h+=noise(x*0.02f,z*0.02f)*10;
    h+=noise(x*0.1f,z*0.1f)*3;
    if(h<1)h=1;
    if(h>80)h=80;
    return h;
}

Chunk* find_chunk(int cx,int cz){
    for(int i=0;i<chunk_count;i++)
        if(chunks[i].cx==cx && chunks[i].cz==cz) return &chunks[i];
    return NULL;
}

Chunk* make_chunk(int cx,int cz){
    if(chunk_count>=MAX_CHUNKS){
        fprintf(stderr,"Too many chunks\n");
        exit(1);
    }
    Chunk* c=&chunks[chunk_count++];
    c->cx=cx;
    c->cz=cz;
    c->count=0;
    c->blocks=malloc(sizeof(Block)*MAX_BLOCKS_PER_CHUNK);
    for(int x=0;x<CHUNK_SIZE;x++)
        for(int y=0;y<CHUNK_HEIGHT;y++)
            for(int z=0;z<CHUNK_SIZE;z++)
                c->map[x][y][z]=MAT_AIR;
    c->dlist=0;
    return c;
}

void add_block_to_map(Chunk* c,int lx,int ly,int lz,Material mat){
    if(lx>=0&&lx<CHUNK_SIZE&&ly>=0&&ly<CHUNK_HEIGHT&&lz>=0&&lz<CHUNK_SIZE&&c->count<MAX_BLOCKS_PER_CHUNK){
        float wx_corner=(c->cx*CHUNK_SIZE+lx)*BLOCK_SIZE;
        float wz_corner=(c->cz*CHUNK_SIZE+lz)*BLOCK_SIZE;
        float wy_corner=ly*BLOCK_SIZE;
        float wx=wx_corner+BLOCK_SIZE/2.0f;
        float wy=wy_corner+BLOCK_SIZE/2.0f;
        float wz=wz_corner+BLOCK_SIZE/2.0f;
        c->blocks[c->count++] = (Block){wx,wy,wz,mat};
        c->map[lx][ly][lz] = mat;
    }
}

#define IS_SOLID(nx,ny,nz) (nx>=0 && nx<CHUNK_SIZE && ny>=0 && ny<CHUNK_HEIGHT && nz>=0 && nz<CHUNK_SIZE && material_table[c->map[nx][ny][nz]].solid)

void draw_face(float x1,float y1,float z1,float x2,float y2,float z2,float x3,float y3,float z3,float x4,float y4,float z4){
    glBegin(GL_QUADS);
    glVertex3f(x1,y1,z1);
    glVertex3f(x2,y2,z2);
    glVertex3f(x3,y3,z3);
    glVertex3f(x4,y4,z4);
    glEnd();
}

void draw_block_optimized(const Block* b,int lx,int ly,int lz,const Chunk* c){
    MaterialInfo m=material_table[b->mat];
    glColor3f(m.r,m.g,m.b);
    float X=b->x,Y=b->y,Z=b->z;
    float S=BLOCK_SIZE;
    if(!IS_SOLID(lx,ly+1,lz)){glNormal3f(0,1,0);draw_face(X-S/2,Y+S/2,Z-S/2,X+S/2,Y+S/2,Z-S/2,X+S/2,Y+S/2,Z+S/2,X-S/2,Y+S/2,Z+S/2);}
    if(!IS_SOLID(lx,ly-1,lz)){glNormal3f(0,-1,0);draw_face(X-S/2,Y-S/2,Z-S/2,X-S/2,Y-S/2,Z+S/2,X+S/2,Y-S/2,Z+S/2,X+S/2,Y-S/2,Z-S/2);}
    if(!IS_SOLID(lx,ly,lz+1)){glNormal3f(0,0,1);draw_face(X-S/2,Y-S/2,Z+S/2,X+S/2,Y-S/2,Z+S/2,X+S/2,Y+S/2,Z+S/2,X-S/2,Y+S/2,Z+S/2);}
    if(!IS_SOLID(lx,ly,lz-1)){glNormal3f(0,0,-1);draw_face(X-S/2,Y-S/2,Z-S/2,X-S/2,Y+S/2,Z-S/2,X+S/2,Y+S/2,Z-S/2,X+S/2,Y-S/2,Z-S/2);}
    if(!IS_SOLID(lx+1,ly,lz)){glNormal3f(1,0,0);draw_face(X+S/2,Y-S/2,Z-S/2,X+S/2,Y+S/2,Z-S/2,X+S/2,Y+S/2,Z+S/2,X+S/2,Y-S/2,Z+S/2);}
    if(!IS_SOLID(lx-1,ly,lz)){glNormal3f(-1,0,0);draw_face(X-S/2,Y-S/2,Z-S/2,X-S/2,Y-S/2,Z+S/2,X-S/2,Y+S/2,Z+S/2,X-S/2,Y+S/2,Z-S/2);}
}

void generate_tree(Chunk* c,int lx,int ly,int lz){
    int h=3+(hash_coords((c->cx*CHUNK_SIZE+lx),(c->cz*CHUNK_SIZE+lz))%3);
    for(int i=0;i<h;i++)
        add_block_to_map(c,lx,ly+i,lz,MAT_WOOD);
    for(int dx=-2;dx<=2;dx++)
        for(int dz=-2;dz<=2;dz++)
            for(int dy=0;dy<=2;dy++){
                if(fabs(dx)+fabs(dz)+fabs(dy)<4)
                    add_block_to_map(c,lx+dx,ly+h+dy,lz+dz,MAT_LEAVES);
            }
}

void generate_flower(Chunk* c,int lx,int ly,int lz){
    int id=hash_coords((c->cx*CHUNK_SIZE+lx),(c->cz*CHUNK_SIZE+lz))%3;
    if(id==0)
        add_block_to_map(c,lx,ly,lz,MAT_FLOWER_RED);
    else if(id==1)
        add_block_to_map(c,lx,ly,lz,MAT_FLOWER_YELLOW);
    else
        add_block_to_map(c,lx,ly,lz,MAT_FLOWER_BLUE);
}

void generate_animal_block(Chunk* c,int lx,int ly,int lz){
    add_block_to_map(c,lx,ly,lz,MAT_ANIMAL);
}

void generate_house_fixed(Chunk* c,int lx,int ly,int lz){
    int w=5+(hash_coords((c->cx*CHUNK_SIZE+lx+1),(c->cz*CHUNK_SIZE+lz))%3);
    int l=5+(hash_coords((c->cx*CHUNK_SIZE+lx+2),(c->cz*CHUNK_SIZE+lz))%3);
    int h=4+(hash_coords((c->cx*CHUNK_SIZE+lx+3),(c->cz*CHUNK_SIZE+lz))%2);
    int base_y=ly;
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            if(gy>base_y) base_y=gy;
        }
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            for(int fy=gy;fy<=base_y;fy++)
                add_block_to_map(c,lx+px,fy,lz+pz,MAT_STONE);
            for(int py=0;py<h;py++){
                int by=base_y+1+py;
                int wall=(px==0||px==w-1||pz==0||pz==l-1||py==0||py==h-1);
                int door=(px==w/2&&pz==0&&py<2);
                int window=((px==0||px==w-1||pz==0||pz==l-1)&&py==2&&(px%2==0||pz%2==0));
                if(wall&&!door&&!window)
                    add_block_to_map(c,lx+px,by,lz+pz,MAT_WOOD);
            }
        }
    for(int level=0;level<3;level++)
        for(int px=-level;px<w+level;px++)
            for(int pz=-level;pz<l+level;pz++)
                add_block_to_map(c,lx+px,base_y+1+h+level,lz+pz,MAT_ROOF);
}

void generate_house_cottage(Chunk* c,int lx,int ly,int lz){
    int w=5+(hash_coords(lx+1,lz)%3);
    int l=5+(hash_coords(lx+2,lz)%3);
    int h=4+(hash_coords(lx+3,lz)%2);
    int base_y=ly;
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            if(gy>base_y) base_y=gy;
        }
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            for(int fy=gy;fy<=base_y;fy++)
                add_block_to_map(c,lx+px,fy,lz+pz,MAT_STONE);
            for(int py=0;py<h;py++){
                int by=base_y+1+py;
                int wall=(px==0||px==w-1||pz==0||pz==l-1||py==0||py==h-1);
                int door=(px==w/2&&pz==0&&py<2);
                int window=((px==0||px==w-1||pz==0||pz==l-1)&&py==2&&(px%2==0||pz%2==0));
                if(wall&&!door&&!window)
                    add_block_to_map(c,lx+px,by,lz+pz,MAT_WOOD);
            }
        }
    for(int level=0;level<3;level++)
        for(int px=-level;px<w+level;px++)
            for(int pz=-level;pz<l+level;pz++)
                add_block_to_map(c,lx+px,base_y+1+h+level,lz+pz,MAT_ROOF);
}

void generate_house_hut(Chunk* c,int lx,int ly,int lz){
    int w=4,l=4,h=3;
    int base_y=ly;
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            if(gy>base_y) base_y=gy;
        }
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            for(int py=0;py<h;py++){
                int wall=(px==0||px==w-1||pz==0||pz==l-1);
                if(wall)
                    add_block_to_map(c,lx+px,base_y+1+py,lz+pz,MAT_BRICK);
            }
        }
    for(int px=-1;px<w+1;px++)
        for(int pz=-1;pz<l+1;pz++)
            add_block_to_map(c,lx+px,base_y+1+h,lz+pz,MAT_STONE);
}

void generate_house_tower(Chunk* c,int lx,int ly,int lz){
    int w=3,l=3,h=8;
    int base_y=ly;
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++){
            float wx=(c->cx*CHUNK_SIZE+lx+px)*BLOCK_SIZE;
            float wz=(c->cz*CHUNK_SIZE+lz+pz)*BLOCK_SIZE;
            int gy=(int)floor(get_ground_height(wx,wz)/BLOCK_SIZE);
            if(gy>base_y) base_y=gy;
        }
    for(int px=0;px<w;px++)
        for(int pz=0;pz<l;pz++)
            for(int py=0;py<h;py++){
                int wall=(px==0||px==w-1||pz==0||pz==l-1);
                if(wall)
                    add_block_to_map(c,lx+px,base_y+1+py,lz+pz,MAT_STONE);
            }
    for(int px=-1;px<w+1;px++)
        for(int pz=-1;pz<l+1;pz++)
            add_block_to_map(c,lx+px,base_y+1+h,lz+pz,MAT_STONE);
}

void generate_house_variant(Chunk* c,int lx,int ly,int lz){
    int choice=hash_coords(lx,lz)%3;
    if(choice==0)
        generate_house_cottage(c,lx,ly,lz);
    else if(choice==1)
        generate_house_hut(c,lx,ly,lz);
    else
        generate_house_tower(c,lx,ly,lz);
}

void generate_village(Chunk* c,int lx,int ly,int lz){
    int num_houses=3+(hash_coords(lx,lz)%4);
    int radius=8;
    for(int i=0;i<num_houses;i++){
        int ox=((int)hash_coords(lx+i,lz+i)% (radius*2))-radius;
        int oz=((int)hash_coords(lz+i,lx+i)% (radius*2))-radius;
        int nx=lx+ox;
        int nz=lz+oz;
        generate_house_variant(c,nx,ly,nz);
    }
}
void reserve_space(int cx, int cz, int lx, int lz, int width, int length, int reserved_map[CHUNK_SIZE][CHUNK_SIZE]){
    int margin = 3; 
    int min_x = lx - margin;
    int max_x = lx + width + margin;
    int min_z = lz - margin;
    int max_z = lz + length + margin;
    if (min_x < 0) min_x = 0;
    if (max_x >= CHUNK_SIZE) max_x = CHUNK_SIZE - 1;
    if (min_z < 0) min_z = 0;
    if (max_z >= CHUNK_SIZE) max_z = CHUNK_SIZE - 1;
    for (int x = min_x; x <= max_x; x++) {
        for (int z = min_z; z <= max_z; z++) {
            reserved_map[x][z] = 1;
        }
    }
}
void generate_chunk(int cx, int cz){
    if(find_chunk(cx,cz)) return;
    Chunk* c=make_chunk(cx,cz);
    int reserved_map[CHUNK_SIZE][CHUNK_SIZE] = {0}; 
    
    for(int x=0;x<CHUNK_SIZE;x++)
        for(int z=0;z<CHUNK_SIZE;z++){
            float wx_sample=(cx*CHUNK_SIZE+x)*BLOCK_SIZE;
            float wz_sample=(cz*CHUNK_SIZE+z)*BLOCK_SIZE;
            float h=get_ground_height(wx_sample,wz_sample);
            int ih=(int)floor(h/BLOCK_SIZE);
            for(int y=0;y<=ih;y++){
                Material m;
                if(y==ih) m=MAT_GRASS;
                else if(y>ih-4) m=MAT_DIRT;
                else m=MAT_STONE;
                add_block_to_map(c,x,y,z,m);
            }
            
            unsigned hv=hash_coords((int)wx_sample,(int)wz_sample);
            int sy=ih+1;
            if(sy<CHUNK_HEIGHT && reserved_map[x][z] == 0){ 
                
                if(hv%300==0) {
                    reserve_space(cx, cz, x, z, 16, 16, reserved_map);
                    generate_village(c,x,sy,z);
                }
                else if(hv%80==0) {
                    reserve_space(cx, cz, x, z, 8, 8, reserved_map); 
                    generate_house_variant(c,x,sy,z);
                } 
            }
            if(sy<CHUNK_HEIGHT && reserved_map[x][z] == 0){ 
                if(hv%40==0)
                    generate_tree(c,x,sy,z);
                else if(hv%20==0)
                    generate_flower(c,x,sy,z);
                else if(hv%100==0)
                    generate_animal_block(c,x,sy,z);
            }
        }
    c->dlist=glGenLists(1);
    glNewList(c->dlist,GL_COMPILE);
    for(int i=0;i<c->count;i++){
        int lx=(int)roundf(c->blocks[i].x/BLOCK_SIZE-0.5f)-c->cx*CHUNK_SIZE;
        int lz=(int)roundf(c->blocks[i].z/BLOCK_SIZE-0.5f)-c->cz*CHUNK_SIZE;
        int ly=(int)roundf(c->blocks[i].y/BLOCK_SIZE-0.5f);
        draw_block_optimized(&c->blocks[i],lx,ly,lz,c);
    }
    glEndList();
    free(c->blocks);
    c->blocks=NULL;
}

void render_text(float x,float y,const char*str){
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,win_w,0,win_h);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_LIGHTING);
    glColor3f(1,1,1);
    glRasterPos2f(x,y);
    for(const char*p=str;*p;p++)
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*p);
    glEnable(GL_LIGHTING);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

void draw_hud(){
    char buf[128];
    int day_seconds=(int)fmod(world_time,DAY_LENGTH);
    int hours=(int)((day_seconds/(float)DAY_LENGTH)*24.0f);
    int minutes=(int)(((day_seconds/(float)DAY_LENGTH)*24.0f-hours)*60.0f);
    snprintf(buf,sizeof(buf),"Time: %02d:%02d",hours,minutes);
    render_text(10,win_h-20,buf);
    snprintf(buf,sizeof(buf),"XYZ: %.1f %.1f %.1f",player_x,player_y,player_z);
    render_text(10,win_h-40,buf);
}

void draw_menu(){
    render_text(win_w/2-30,win_h/2+20,"PAUSED");
    render_text(win_w/2-60,win_h/2,"R - Resume");
    render_text(win_w/2-60,win_h/2-20,"Q - Quit");
}

void ensure_chunks_fixed(float px,float pz){
    int pcx=(int)floor(px/(CHUNK_SIZE*BLOCK_SIZE)), pcz=(int)floor(pz/(CHUNK_SIZE*BLOCK_SIZE));
    for(int dx=-WORLD_RADIUS;dx<=WORLD_RADIUS;dx++)
        for(int dz=-WORLD_RADIUS;dz<=WORLD_RADIUS;dz++)
            if(!find_chunk(pcx+dx,pcz+dz))
                generate_chunk(pcx+dx,pcz+dz);
}

void draw_world(){
    for(int i=0;i<chunk_count;i++)
        glCallList(chunks[i].dlist);
}

void update_lighting(float dt){
    world_time+=dt;
    if(world_time>DAY_LENGTH) world_time-=DAY_LENGTH;
    float time_factor=sinf(world_time/DAY_LENGTH*2.0f*M_PI);
    time_factor=(time_factor*0.5f)+0.5f;
    float min_light=0.2f, max_light=1.0f;
    float light_intensity=min_light+(max_light-min_light)*time_factor;
    float sky_day[]={0.5f,0.7f,1.0f};
    float sky_night[]={0.05f,0.05f,0.15f};
    float sky_r=sky_night[0]+(sky_day[0]-sky_night[0])*time_factor;
    float sky_g=sky_night[1]+(sky_day[1]-sky_night[1])*time_factor;
    float sky_b=sky_night[2]+(sky_day[2]-sky_night[2])*time_factor;
    glClearColor(sky_r,sky_g,sky_b,1.0f);
    GLfloat ambient[]={light_intensity*0.7f,light_intensity*0.7f,light_intensity*0.7f,1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambient);
    GLfloat diffuse[]={light_intensity,light_intensity,light_intensity,1.0f};
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);
    float angle=world_time/DAY_LENGTH*2.0f*M_PI;
    GLfloat light_pos[]={sinf(angle)*50.0f,cosf(angle)*50.0f,0.0f,0.0f};
    glLightfv(GL_LIGHT0,GL_POSITION,light_pos);
}

int block_at(int bx,int by,int bz){
    int cx=bx/CHUNK_SIZE;
    if(bx<0 && bx%CHUNK_SIZE) cx--;
    int cz=bz/CHUNK_SIZE;
    if(bz<0 && bz%CHUNK_SIZE) cz--;
    Chunk* c=find_chunk(cx,cz);
    if(!c) return MAT_AIR;
    int lx=(bx%CHUNK_SIZE+CHUNK_SIZE)%CHUNK_SIZE;
    int lz=(bz%CHUNK_SIZE+CHUNK_SIZE)%CHUNK_SIZE;
    if(by<0||by>=CHUNK_HEIGHT) return MAT_AIR;
    return c->map[lx][by][lz];
}

int collide_player(float nx,float ny,float nz){
    int minx=(int)floor((nx-PLAYER_RADIUS)/BLOCK_SIZE), maxx=(int)floor((nx+PLAYER_RADIUS)/BLOCK_SIZE);
    int miny=(int)floor(ny/BLOCK_SIZE), maxy=(int)floor((ny+PLAYER_HEIGHT)/BLOCK_SIZE);
    int minz=(int)floor((nz-PLAYER_RADIUS)/BLOCK_SIZE), maxz=(int)floor((nz+PLAYER_RADIUS)/BLOCK_SIZE);
    for(int x=minx;x<=maxx;x++)
        for(int y=miny;y<=maxy;y++)
            for(int z=minz;z<=maxz;z++){
                int mat=block_at(x,y,z);
                if(mat!=MAT_AIR && material_table[mat].solid) return 1;
            }
    return 0;
}

void move_player_axis(float dx,float dy,float dz){
    if(dx!=0){
        if(!collide_player(player_x+dx,player_y,player_z))
            player_x+=dx;
    }
    if(dz!=0){
        if(!collide_player(player_x,player_y,player_z+dz))
            player_z+=dz;
    }
    if(dy!=0){
        if(!collide_player(player_x,player_y+dy,player_z))
            player_y+=dy;
    }
}

void update_player(float dt){
    if(!paused) update_lighting(dt);
    float fx=0,fz=0;
    float speed=is_flying?FLY_SPEED:MOVE_SPEED;
    if(key_w) fx+=sinf(yaw*M_PI/180.0f), fz+=cosf(yaw*M_PI/180.0f);
    if(key_s) fx-=sinf(yaw*M_PI/180.0f), fz-=cosf(yaw*M_PI/180.0f);
    if(key_a) fx-=cosf(yaw*M_PI/180.0f), fz+=sinf(yaw*M_PI/180.0f);
    if(key_d) fx+=cosf(yaw*M_PI/180.0f), fz-=sinf(yaw*M_PI/180.0f);
    float len2=fx*fx+fz*fz;
    if(len2>0){
        float inv_len=q_rsqrt(len2);
        fx*=inv_len;
        fz*=inv_len;
    }
    float dx=fx*speed*dt;
    float dz=fz*speed*dt;
    move_player_axis(dx,0,dz);
    if(is_flying){
        if(key_space) move_player_axis(0,FLY_SPEED*dt,0);
        if(key_c) move_player_axis(0,-FLY_SPEED*dt,0);
        vel_y=0;
        on_ground=0;
    }
    else {
        vel_y+=GRAVITY*dt;
        if(on_ground && key_space){
            vel_y=JUMP_SPEED;
            on_ground=0;
        }
        float try_y=player_y+vel_y*dt;
        if(!collide_player(player_x,try_y,player_z))
            player_y=try_y;
        else {
            if(vel_y<0){
                int foot_block=(int)floor(player_y/BLOCK_SIZE);
                while(foot_block>=-1){
                    if(block_at((int)floor(player_x/BLOCK_SIZE),foot_block,(int)floor(player_z/BLOCK_SIZE))!=MAT_AIR){
                        player_y=(foot_block+1)*BLOCK_SIZE;
                        vel_y=0;
                        on_ground=1;
                        break;
                    }
                    foot_block--;
                }
            }
            else vel_y=0;
        }
        float groundh=get_ground_height(player_x,player_z);
        float ground_floor=floor(groundh/BLOCK_SIZE)*BLOCK_SIZE;
        if(player_y<ground_floor+0.05f){
            player_y=ground_floor+0.05f;
            vel_y=0;
            on_ground=1;
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
    gluLookAt(player_x,player_y+PLAYER_HEIGHT*0.5f,player_z,
              player_x+dx,player_y+PLAYER_HEIGHT*0.5f+dy,player_z+dz,0,1,0);
    draw_world();
    draw_hud();
    if(paused) draw_menu();
    glutSwapBuffers();
    glutPostRedisplay();
}

void reshape(int w,int h){
    win_w=w;
    win_h=h;
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(65,(float)w/h,0.1,1000);
    glMatrixMode(GL_MODELVIEW);
    glutWarpPointer(win_w/2, win_h/2);
}

void key_down(unsigned char k,int x,int y){
    k=tolower(k);
    if(paused){
        if(k=='r') paused=0;
        if(k=='q') exit(0);
        return;
    }
    if(k==27){ paused=1; return;}
    if(k=='w') key_w=1;
    if(k=='s') key_s=1;
    if(k=='a') key_a=1;
    if(k=='d') key_d=1;
    if(k==' ') key_space=1;
    if(k=='c') key_c=1;
    if(k=='f') is_flying=!is_flying;
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

int find_village_chunk_coord(int start_coord){
    for(int i=0;i<500;i++){
        int cx=start_coord+((hash(i)%30)-15);
        int cz=start_coord+((hash(i+1)%30)-15);
        if(hash_coords(cx*CHUNK_SIZE, cz*CHUNK_SIZE)%300 == 0)
            return cx;
    }
    return start_coord;
}

int main(int argc,char**argv){
    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    int screen_w=glutGet(GLUT_SCREEN_WIDTH);
    int screen_h=glutGet(GLUT_SCREEN_HEIGHT);
    glutInitWindowSize(screen_w,screen_h);
    glutInitWindowPosition(0,0);
    glutCreateWindow("minceraft");
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_SMOOTH);
    update_lighting(0.0f);
    glutSetCursor(GLUT_CURSOR_NONE);
    glutWarpPointer(screen_w/2,screen_h/2);

    int initial_cx = find_village_chunk_coord(0);
    int initial_cz = find_village_chunk_coord(0);

    player_x = initial_cx * CHUNK_SIZE * BLOCK_SIZE + CHUNK_SIZE * BLOCK_SIZE / 2.0f;
    player_z = initial_cz * CHUNK_SIZE * BLOCK_SIZE + CHUNK_SIZE * BLOCK_SIZE / 2.0f;

    ensure_chunks_fixed(player_x,player_z);
    float initial_height=get_ground_height(player_x,player_z);
    int initial_block_y=(int)floor(initial_height/BLOCK_SIZE);
    player_y=initial_block_y*BLOCK_SIZE+0.1f;

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key_down);
    glutKeyboardUpFunc(key_up);
    glutPassiveMotionFunc(mouse_motion);
    glutMainLoop();
    for(int i=0;i<chunk_count;i++)
        if(chunks[i].dlist) glDeleteLists(chunks[i].dlist,1);
    return 0;
}