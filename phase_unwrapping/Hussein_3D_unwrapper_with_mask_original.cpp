//This program is written by Hussein Abdul-Rahman and Munther Gdeisat to program the three-dimensional phase unwrapper
//entitled "Fast three-dimensional phase-unwrapping algorithm based on sorting by 
//reliability following a noncontinuous path"
//by  Hussein Abdul-Rahman, Munther A. Gdeisat, David R. Burton, and Michael J. Lalor, 
//published in the Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.
//This program was written on 29th August 2007
//The wrapped phase volume is floating point data type. Also, the unwrapped phase volume is foloating point
//read the data from the file frame by frame
//The mask is byte data type. 
//When the mask is 255 this means that the voxel is valid 
//When the mask is 0 this means that the voxel is invalid (noisy or corrupted voxel)

#include <malloc.h>
// #include "stdafx.h" 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>

static float PI = 3.141592654;
static float TWOPI = 6.283185307;
int No_of_edges = 0;

//VOXELM information
struct VOXELM
{
	//int x;					//x coordinate of the voxel
    //int y;					//y coordinate
	//int z;					//z coordinate
    int increment;				//No. of 2*pi to add to the voxel to unwrap it
    int number_of_voxels_in_group;	//No. of voxel in the voxel group
    float value;				//value of the voxel
	float reliability;
	unsigned char input_mask;			//0 voxel is masked. 255 voxel is not masked
	unsigned char extended_mask;			//0 voxel is masked. 255 voxel is not masked
    int group;					//group No.
    int new_group;
    struct VOXELM *head;		//pointer to the first voxel in the group in the linked list
    struct VOXELM *last;		//pointer to the last voxel in the group
    struct VOXELM *next;		//pointer to the next voxel in the group
};


//the EDGE is the line that connects two voxels.
//if we have S voxels, then we have S horizental edges and S vertical edges
struct EDGE
{    
	float reliab;			//reliabilty of the edge and it depends on the two voxels
	VOXELM *pointer_1;		//pointer to the first voxel
    VOXELM *pointer_2;		//pointer to the second voxel
    int increment;			//No. of 2*pi to add to one of the voxels to unwrap it with respect to the second 
}; 

void read_data(char *inputfile,float *Data, int length)
{
	printf("Reading the Wrapped Values form Binary File.............>");
	FILE *ifptr;
	ifptr = fopen(inputfile,"rb");
	if(ifptr == NULL) printf("Error opening the file\n");
	fread(Data,sizeof(float),length,ifptr);
	fclose(ifptr);
	printf(" Done.\n");
}

void write_data(char *outputfile,float *Data,int length)
{
	printf("Writing the Unwrapped Values to Binary File.............>");
	FILE *ifptr;
	ifptr = fopen(outputfile,"wb");
	if(ifptr == NULL) printf("Error opening the file\n");
	fwrite(Data,sizeof(float),length,ifptr);
	fclose(ifptr);
	printf(" Done.\n");
}

void read_mask(char *inputfile,unsigned char *Data, int length)
{
	printf("Reading the mask form Binary File.............>");
	FILE *ifptr;
	ifptr = fopen(inputfile,"rb");
	if(ifptr == NULL) printf("Error opening the file\n");
	fread(Data,sizeof(char),length,ifptr);
	fclose(ifptr);
	printf(" Done.\n");
}

//another version of Mixtogether but this function should only be use with the sort program
void  Mix(EDGE *Pointer1, int *index1, int *index2, int size)
{
	int counter1 = 0;
	int counter2 = 0;
	int *TemporalPointer = index1;

	int *Result = (int *) calloc(size * 2, sizeof(int));
	int *Follower = Result;

	while ((counter1 < size) && (counter2 < size))
	{
		if ((Pointer1[*(index1 + counter1)].reliab <= Pointer1[*(index2 + counter2)].reliab))
		{
			*Follower = *(index1 + counter1);
			Follower++;
			counter1++;
		} 
		else
        {
			*Follower = *(index2 + counter2);
			Follower++;
			counter2++;
        }
	}//while

	if (counter1 == size)
	{
		memcpy(Follower, (index2 + counter2), sizeof(int)*(size-counter2));
	} 
	else
	{
		memcpy(Follower, (index1 + counter1), sizeof(int)*(size-counter1));
	}

	Follower = Result;
	index1 = TemporalPointer;

	int i;
	for (i=0; i < 2 * size; i++)
	{
		*index1 = *Follower;
		index1++;
		Follower++;
	}

	free(Result);
}

//this is may be the fastest sort program; 
//see the explination in quickSort function below
void  sort(EDGE *Pointer, int *index, int size)
{
	if (size == 2)
	{
		if ((Pointer[*index].reliab) > (Pointer[*(index+1)].reliab))
		{
			int Temp;
			Temp = *index;
			*index = *(index+1);
			*(index+1) = Temp;
		}
	} 
	else if (size > 2)
    {
		sort(Pointer, index, size/2);
		sort(Pointer, (index + (size/2)), size/2);
		Mix(Pointer, index, (index + (size/2)), size/2);
    }
}

//this function tries to implement a nice idea explained below
//we need to sort edge array. Each edge element conisists of 16 bytes.
//In normal sort program we compare two elements in the array and exchange
//their place under some conditions to do the sorting. It is very probable
// that an edge element may change its place hundred of times which makes 
//the sorting a very time consuming operation. The idea in this function 
//is to give each edge element an index and move the index not the edge
//element. The edge need 4 bytes which makes the sorting operation faster.
// After finishingthe sorting of the indexes, we know the position of each index. 
//So we know how to sort edges
void  quick_sort(EDGE *Pointer, int size)
{
	int *index = (int *) calloc(size, sizeof(int));
	int i;

	for (i=0; i<size; ++i)
		index[i] = i;

	sort(Pointer, index, size);

	EDGE * a = (EDGE *) calloc(size, sizeof(EDGE));
	for (i=0; i<size; ++i)
		a[i] = Pointer[*(index + i)];

	memcpy(Pointer, a, size*sizeof(EDGE));

	free(index);
	free(a);
}

//---------------start quicker_sort algorithm --------------------------------
#define swap(x,y) {EDGE t; t=x; x=y; y=t;}
#define order(x,y) if (x.reliab > y.reliab) swap(x,y)
#define o2(x,y) order(x,y)
#define o3(x,y,z) o2(x,y); o2(x,z); o2(y,z)

typedef enum {yes, no} yes_no;

yes_no find_pivot(EDGE *left, EDGE *right, float *pivot_ptr)
{
	EDGE a, b, c, *p;

	a = *left;
	b = *(left + (right - left) /2 );
	c = *right;
	o3(a,b,c);

	if (a.reliab < b.reliab)
	{
		*pivot_ptr = b.reliab;
		return yes;
	}

	if (b.reliab < c.reliab)
	{
		*pivot_ptr = c.reliab;
		return yes;
	}

	for (p = left + 1; p <= right; ++p)
	{
		if (p->reliab != left->reliab)
		{
			*pivot_ptr = (p->reliab < left->reliab) ? left->reliab : p->reliab;
			return yes;
		}
		return no;
	}
}

EDGE *partition(EDGE *left, EDGE *right, float pivot)
{
	while (left <= right)
	{
		while (left->reliab < pivot)
			++left;
		while (right->reliab >= pivot)
			--right;
		if (left < right)
		{
			swap (*left, *right);
			++left;
			--right;
		}
	}
	return left;
}

void quicker_sort(EDGE *left, EDGE *right)
{
	EDGE *p;
	float pivot;

	if (find_pivot(left, right, &pivot) == yes)
	{
		p = partition(left, right, pivot);
		quicker_sort(left, p - 1);
		quicker_sort(p, right);
	}
}

//--------------end quicker_sort algorithm -----------------------------------

//--------------------start initialse voxels ----------------------------------
//initialse voxels. See the explination of the voxel class above.
//initially every voxel is a gorup by its self
void  initialiseVOXELs(float *WrappedVolume, unsigned char *input_mask, unsigned char *extended_mask, VOXELM *voxel, int volume_width, int volume_height, int volume_depth)
{
	VOXELM *voxel_pointer = voxel;
	float *wrapped_volume_pointer = WrappedVolume;
	unsigned char *input_mask_pointer = input_mask;
	unsigned char *extended_mask_pointer = extended_mask;
	int n, i, j;

    for (n=0; n < volume_depth; n++)
	{
		for (i=0; i < volume_height; i++)
        {
			for (j=0; j < volume_width; j++)
			{
				//voxel_pointer->x = j;
  				//voxel_pointer->y = i;
				//voxel_pointer->z = n;
				voxel_pointer->increment = 0;
				voxel_pointer->number_of_voxels_in_group = 1;		
  				voxel_pointer->value = *wrapped_volume_pointer;
				voxel_pointer->reliability = 9999999 + rand();
				voxel_pointer->input_mask = *input_mask_pointer;
				voxel_pointer->extended_mask = *extended_mask_pointer;
				voxel_pointer->head = voxel_pointer;
  				voxel_pointer->last = voxel_pointer;
				voxel_pointer->next = NULL;			
				voxel_pointer->new_group = 0;
				voxel_pointer->group = -1;
				voxel_pointer++;
				wrapped_volume_pointer++;
				input_mask_pointer++;
				extended_mask_pointer++;
			}
        }
	}
}
//-------------------end initialise voxels -----------

//gamma finction in the paper
float wrap(float voxel_value)
{
	float wrapped_voxel_value;
	if (voxel_value > PI)	wrapped_voxel_value = voxel_value - TWOPI;
	else if (voxel_value < -PI)	wrapped_voxel_value = voxel_value + TWOPI;
	else wrapped_voxel_value = voxel_value;
	return wrapped_voxel_value;
}

// voxelL_value is the left voxel,	voxelR_value is the right voxel
int find_wrap(float voxelL_value, float voxelR_value)
{
	float difference; 
	int wrap_value;
	difference = voxelL_value - voxelR_value;

	if (difference > PI)	wrap_value = -1;
	else if (difference < -PI)	wrap_value = 1;
	else wrap_value = 0;

	return wrap_value;
} 

void extend_mask(unsigned char *input_mask, unsigned char *extended_mask, int volume_width, int volume_height, int volume_depth)
{
	int n, i, j;
	int vw = volume_width, vh = volume_height, vd = volume_depth;
	int fs  = volume_width * volume_height;	//frame size
	int frame_size = volume_width * volume_height; 
	int volume_size = volume_width * volume_height * volume_depth;	//volume size
	int vs = volume_size;
	unsigned char *IMP = input_mask    + frame_size + volume_width + 1;	//input mask pointer
	unsigned char *EMP = extended_mask + frame_size + volume_width + 1;	//extended mask pointer

	//extend the mask for the volume except borders
	for (n=1; n < volume_depth - 1; n++)
	{
		for (i=1; i < volume_height - 1; i++)
        {
			for (j=1; j < volume_width - 1; j++)
			{
				if( (*IMP) == 255				&& (*(IMP - 1) == 255)				&& (*(IMP + 1) == 255) && 
					(*(IMP + vw) == 255)		&& (*(IMP + vw - 1) == 255)			&& (*(IMP + vw + 1) == 255) &&
					(*(IMP - vw) == 255)		&& (*(IMP - vw - 1) == 255)			&& (*(IMP - vw + 1) == 255) &&
					(*(IMP + fs) == 255)		&& (*(IMP + fs - 1) == 255)			&& (*(IMP + fs + 1) == 255) &&
					(*(IMP + fs - vw) == 255)	&& (*(IMP + fs - vw - 1) == 255)	&& (*(IMP + fs - vw + 1) == 255) &&
					(*(IMP + fs + vw) == 255)	&& (*(IMP + fs + vw - 1) == 255)	&& (*(IMP + fs + vw + 1) == 255) &&
					(*(IMP - fs) == 255)		&& (*(IMP - fs - 1) == 255)			&& (*(IMP - fs + 1) == 255) &&
					(*(IMP - fs - vw) == 255)	&& (*(IMP - fs - vw - 1) == 255)	&& (*(IMP - fs - vw + 1) == 255) &&
					(*(IMP - fs + vw) == 255)	&& (*(IMP - fs + vw - 1) == 255)	&& (*(IMP - fs + vw + 1) == 255))
				{		
					*EMP = 255;
				}
				++EMP;
				++IMP;
			}
			EMP += 2;
			IMP += 2;
		}
		EMP += 2 * volume_width;
		IMP += 2 * volume_width;
	}		
}

void calculate_reliability(float *wrappedVolume, VOXELM *voxel, int volume_width, int volume_height, int volume_depth)
{ 
	int frame_size  = volume_width * volume_height;
	int volume_size = volume_width * volume_height * volume_width;
	VOXELM *voxel_pointer;
	float H, V, N, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10;
	float *WVP;
	int n, i, j;
	
	WVP = wrappedVolume + frame_size + volume_width + 1;
	voxel_pointer = voxel + frame_size + volume_width + 1;
	for (n=1; n < volume_depth - 1; n++)
	{
		for (i=1; i < volume_height - 1; i++)
        {
			for (j=1; j < volume_width - 1; j++)
			{
				if (voxel_pointer->extended_mask == 255)
				{ 
					H  = wrap(*(WVP - 1) - *WVP) - wrap(*WVP - *(WVP + 1));
					V  = wrap(*(WVP - volume_width) - *WVP) - wrap(*WVP - *(WVP + volume_width));
					N  = wrap(*(WVP - frame_size) - *WVP) - wrap(*WVP - *(WVP + frame_size));
					D1 = wrap(*(WVP - volume_width - 1) - *WVP) - wrap(*WVP - *(WVP + volume_width + 1));
					D2 = wrap(*(WVP - volume_width + 1) - *WVP) - wrap(*WVP - *(WVP + volume_width - 1));
					D3 = wrap(*(WVP - frame_size - volume_width - 1) - *WVP) - wrap(*WVP - *(WVP + frame_size + volume_width + 1));
					D4 = wrap(*(WVP - frame_size - volume_width) - *WVP) - wrap(*WVP - *(WVP + frame_size + volume_width));
					D5 = wrap(*(WVP - frame_size - volume_width + 1) - *WVP) - wrap(*WVP - *(WVP + frame_size + volume_width - 1));
					D6 = wrap(*(WVP - frame_size - 1) - *WVP) - wrap(*WVP - *(WVP + frame_size + 1));
					D7 = wrap(*(WVP - frame_size + 1) - *WVP) - wrap(*WVP - *(WVP + frame_size - 1));
					D8 = wrap(*(WVP - frame_size + volume_width - 1) - *WVP) - wrap(*WVP - *(WVP + frame_size - volume_width + 1));
					D9 = wrap(*(WVP - frame_size + volume_width) - *WVP) - wrap(*WVP - *(WVP + frame_size - volume_width));
					D10 = wrap(*(WVP - frame_size + volume_width + 1) - *WVP) - wrap(*WVP - *(WVP + frame_size - volume_width - 1));
					voxel_pointer->reliability = H*H + V*V + N*N + D1*D1 + D2*D2  + D3*D3 + D4*D4  + D5*D5 + D6*D6  
						+ D7*D7 + D8*D8 + D9*D9 + D10*D10;
				}
				voxel_pointer++;
				WVP++;
			}
			voxel_pointer += 2;
			WVP += 2;
		}
		voxel_pointer += 2 * volume_width;
		WVP += 2 * volume_width;
	}
}

//calculate the reliability of the horizental edges of the volume
//it is calculated by adding the reliability of voxel and the relibility of 
//its right neighbour
//edge is calculated between a voxel and its next neighbour
void  horizentalEDGEs(VOXELM *voxel, EDGE *edge, int volume_width, int volume_height, int volume_depth)
{
	int n, i, j;
	EDGE *edge_pointer = edge;
	VOXELM *voxel_pointer = voxel;
	
	for (n=0; n < volume_depth; n++)
	{
		for (i = 0; i < volume_height; i++)
		{
			for (j = 0; j < volume_width - 1; j++) 
			{
				if (voxel_pointer->input_mask == 255 && (voxel_pointer + 1)->input_mask == 255 )
				{
					edge_pointer->pointer_1 = voxel_pointer;
					edge_pointer->pointer_2 = (voxel_pointer+1);
					edge_pointer->reliab = voxel_pointer->reliability + (voxel_pointer + 1)->reliability;
					edge_pointer->increment = find_wrap(voxel_pointer->value, (voxel_pointer + 1)->value);
					edge_pointer++;
					No_of_edges++;
				}
				voxel_pointer++;
			}
			voxel_pointer++;
		}
	}
}

void  verticalEDGEs(VOXELM *voxel, EDGE *edge, int volume_width, int volume_height, int volume_depth)
{
	int n, i, j;	
	VOXELM *voxel_pointer = voxel;
	EDGE *edge_pointer = edge + No_of_edges; 

	for (n=0; n < volume_depth; n++)
	{
		for (i=0; i<volume_height - 1; i++)
		{
			for (j=0; j < volume_width; j++) 
			{
				if (voxel_pointer->input_mask == 255 && (voxel_pointer + volume_width)->input_mask == 255 )
				{
					edge_pointer->pointer_1 = voxel_pointer;
					edge_pointer->pointer_2 = (voxel_pointer + volume_width);
					edge_pointer->reliab = voxel_pointer->reliability + (voxel_pointer + volume_width)->reliability;
					edge_pointer->increment = find_wrap(voxel_pointer->value, (voxel_pointer + volume_width)->value);
					edge_pointer++;
					No_of_edges++;
				}
				voxel_pointer++;
			}
		}
		voxel_pointer += volume_width;
	} 
}

void  normalEDGEs(VOXELM *voxel, EDGE *edge, int volume_width, int volume_height, int volume_depth)
{
	int n, i, j;	
	int frame_size = volume_width * volume_height;
	VOXELM *voxel_pointer = voxel;
	EDGE *edge_pointer = edge + No_of_edges; 

	for (n=0; n < volume_depth - 1; n++)
	{
		for (i=0; i<volume_height; i++)
		{
			for (j=0; j < volume_width; j++) 
			{
				if (voxel_pointer->input_mask == 255 && (voxel_pointer + frame_size)->input_mask == 255 )
				{
					edge_pointer->pointer_1 = voxel_pointer;
					edge_pointer->pointer_2 = (voxel_pointer + frame_size);
					edge_pointer->reliab = voxel_pointer->reliability + (voxel_pointer + frame_size)->reliability;
					edge_pointer->increment = find_wrap(voxel_pointer->value, (voxel_pointer + frame_size)->value);
					edge_pointer++;
					No_of_edges++;
				}
				voxel_pointer++;
			}
		}
	} 
}

//gather the voxels of the volume into groups 
void  gatherVOXELs(EDGE *edge)
{
	int k;
	VOXELM *VOXEL1;   
	VOXELM *VOXEL2;
	VOXELM *group1;
	VOXELM *group2;
	EDGE *pointer_edge = edge;
	int incremento;

	for (k = 0; k < No_of_edges; k++)
	{
		VOXEL1 = pointer_edge->pointer_1;
		VOXEL2 = pointer_edge->pointer_2;

		//VOXELM 1 and VOXELM 2 belong to different groups
		//initially each voxel is a group by it self and one voxel can construct a group
		//no else or else if to this if
		if (VOXEL2->head != VOXEL1->head)
		{
			//VOXELM 2 is alone in its group
			//merge this voxel with VOXELM 1 group and find the number of 2 pi to add 
			//to or subtract to unwrap it
			if ((VOXEL2->next == NULL) && (VOXEL2->head == VOXEL2))
			{
				VOXEL1->head->last->next = VOXEL2;
				VOXEL1->head->last = VOXEL2;
				(VOXEL1->head->number_of_voxels_in_group)++;
				VOXEL2->head=VOXEL1->head;
				VOXEL2->increment = VOXEL1->increment-pointer_edge->increment;
			}

			//VOXELM 1 is alone in its group
			//merge this voxel with VOXELM 2 group and find the number of 2 pi to add 
			//to or subtract to unwrap it
			else if ((VOXEL1->next == NULL) && (VOXEL1->head == VOXEL1))
			{
				VOXEL2->head->last->next = VOXEL1;
				VOXEL2->head->last = VOXEL1;
				(VOXEL2->head->number_of_voxels_in_group)++;
				VOXEL1->head = VOXEL2->head;
				VOXEL1->increment = VOXEL2->increment+pointer_edge->increment;
			} 

			//VOXELM 1 and VOXELM 2 both have groups
			else
            {
				group1 = VOXEL1->head;
                group2 = VOXEL2->head;
				//the no. of voxels in VOXELM 1 group is large than the no. of voxels
				//in VOXELM 2 group.   Merge VOXELM 2 group to VOXELM 1 group
				//and find the number of wraps between VOXELM 2 group and VOXELM 1 group
				//to unwrap VOXELM 2 group with respect to VOXELM 1 group.
				//the no. of wraps will be added to VOXELM 2 grop in the future
				if (group1->number_of_voxels_in_group > group2->number_of_voxels_in_group)
				{
					//merge VOXELM 2 with VOXELM 1 group
					group1->last->next = group2;
					group1->last = group2->last;
					group1->number_of_voxels_in_group = group1->number_of_voxels_in_group + group2->number_of_voxels_in_group;
					incremento = VOXEL1->increment-pointer_edge->increment - VOXEL2->increment;
					//merge the other voxels in VOXELM 2 group to VOXELM 1 group
					while (group2 != NULL)
					{
						group2->head = group1;
						group2->increment += incremento;
						group2 = group2->next;
					}
				} 

				//the no. of voxels in VOXELM 2 group is large than the no. of voxels
				//in VOXELM 1 group.   Merge VOXELM 1 group to VOXELM 2 group
				//and find the number of wraps between VOXELM 2 group and VOXELM 1 group
				//to unwrap VOXELM 1 group with respect to VOXELM 2 group.
				//the no. of wraps will be added to VOXELM 1 grop in the future
				else
                {
					//merge VOXELM 1 with VOXELM 2 group
					group2->last->next = group1;
					group2->last = group1->last;
					group2->number_of_voxels_in_group = group2->number_of_voxels_in_group + group1->number_of_voxels_in_group;
					incremento = VOXEL2->increment + pointer_edge->increment - VOXEL1->increment;
					//merge the other voxels in VOXELM 2 group to VOXELM 1 group
					while (group1 != NULL)
					{
						group1->head = group2;
						group1->increment += incremento;
						group1 = group1->next;
					} // while

                } // else
            } //else
        } //if
        pointer_edge++;
	}
}

//unwrap the volume 
void  unwrapVolume(VOXELM *voxel, int volume_width, int volume_height, int volume_depth)
{
	int i;
	int volume_size = volume_width * volume_height * volume_depth;
	VOXELM *voxel_pointer=voxel;

	for (i = 0; i < volume_size; i++)
	{
		voxel_pointer->value += TWOPI * (float)(voxel_pointer->increment);
        voxel_pointer++;
    }
}

//set the masked voxels (mask = 0) to the minimum of the unwrapper phase
void  maskVolume(VOXELM *voxel, unsigned char *input_mask, int volume_width, int volume_height, int volume_depth)
{
	int volume_width_plus_one  = volume_width + 1;
	int volume_height_plus_one  = volume_height + 1;
	int volume_width_minus_one = volume_width - 1;
	int volume_height_minus_one = volume_height - 1;

	VOXELM *pointer_voxel = voxel;
	unsigned char *IMP = input_mask;	//input mask pointer
	float min=99999999.;
	int i, j;
	int volume_size = volume_width * volume_height * volume_depth;

	//find the minimum of the unwrapped phase
	for (i = 0; i < volume_size; i++)
	{
		if ((pointer_voxel->value < min) && (*IMP == 255)) 
			min = pointer_voxel->value;

		pointer_voxel++;
		IMP++;
	}

	pointer_voxel = voxel;
	IMP = input_mask;	

	//set the masked voxels to minimum
	for (i = 0; i < volume_size; i++)
	{
		if ((*IMP) == 0)
		{
			pointer_voxel->value = min;
		}
		pointer_voxel++;
		IMP++;
	}
}

//the input to this unwrapper is an array that contains the wrapped phase map. 
//copy the volume on the buffer passed to this unwrapper to over write the unwrapped 
//phase map on the buffer of the wrapped phase map.
void  returnVolume(VOXELM *voxel, float *unwrappedVolume, int volume_width, int volume_height, int volume_depth)
{
	int i;
	int volume_size = volume_width * volume_height * volume_depth;
    float *unwrappedVolume_pointer = unwrappedVolume;
    VOXELM *voxel_pointer = voxel;

    for (i=0; i < volume_size; i++) 
	{
        *unwrappedVolume_pointer = voxel_pointer->value;
        voxel_pointer++;
		unwrappedVolume_pointer++;
	}
}

//the main function of the unwrapper
int main(int argc, char* argv[])
{  
	float *WrappedVolume, *UnwrappedVolume;
	unsigned char *input_mask, *extended_mask;
	int volume_width = std::stoi(argv[4]);
	int volume_height = std::stoi(argv[5]);
	int volume_depth = std::stoi(argv[6]);
	// printf("%d","%d","%d",volume_width,volume_height,volume_depth);
	int volume_size = volume_height * volume_width * volume_depth;
	int No_of_Edges_initially = 3 * volume_width * volume_height * volume_depth;

	WrappedVolume = (float *) calloc(volume_size, sizeof(float));
	read_data(argv[1], WrappedVolume, volume_size);

	input_mask = (unsigned char *) calloc(volume_size, sizeof(unsigned char));
	extended_mask = (unsigned char *) calloc(volume_size, sizeof(unsigned char));
	read_mask(argv[2], input_mask, volume_size);

	UnwrappedVolume = (float *) calloc(volume_size, sizeof(float));
	VOXELM *voxel = (VOXELM *) calloc(volume_size, sizeof(VOXELM));
	EDGE *edge = (EDGE *) calloc(No_of_Edges_initially, sizeof(EDGE));;

	extend_mask(input_mask, extended_mask, volume_width, volume_height, volume_depth);
	
	initialiseVOXELs(WrappedVolume, input_mask, extended_mask, voxel, volume_width, volume_height, volume_depth);

	calculate_reliability(WrappedVolume, voxel, volume_width, volume_height, volume_depth);

	horizentalEDGEs(voxel, edge, volume_width, volume_height, volume_depth);

	verticalEDGEs(voxel, edge, volume_width, volume_height, volume_depth);

	normalEDGEs(voxel, edge, volume_width, volume_height, volume_depth);

	//sort the EDGEs depending on their reiability. The VOXELs with higher relibility (small value) first
	//run only one of the two functions (quick_sort() or quicker_sort() )
	//if your code stuck because of the quicker_sort() function, then use the quick_sort() function
	quick_sort(edge, No_of_edges);
	//quicker_sort(edge, edge + No_of_edges - 1);

	//gather VOXELs into groups
	gatherVOXELs(edge);

	unwrapVolume(voxel, volume_width, volume_height, volume_depth);

	maskVolume(voxel, input_mask, volume_width, volume_height, volume_depth);

	//copy the volume from VOXELM structure to the unwrapped phase array passed to this function
	returnVolume(voxel, UnwrappedVolume, volume_width, volume_height, volume_depth);

	free(edge);
	free(voxel);
	
	write_data(argv[3], UnwrappedVolume, volume_size);
	free(UnwrappedVolume);
	free(WrappedVolume);
	free(input_mask);
	// getchar();
	return 1;
}
