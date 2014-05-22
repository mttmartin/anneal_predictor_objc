//  Created by Matthew Martin on 5/17/14.
//  Copyright (c) 2014 Matthew Martin. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <stdio.h>



char *strrev(const char *string)
{
    char *revstring = calloc(sizeof(char),strlen(string)+1);
    
    unsigned long n = strlen(string);
    
    // go to end of given string
    const char *tmp = string + n - 1;
    
    
    for (int i = 0; i <= n; i++)
    {
        *revstring = *tmp;
        revstring++;
        tmp--;
    }
    
    // put revstring back to initial position
    revstring = revstring - (n) - 1;
    return revstring;
}

@interface DNASequence: NSObject
{
    unsigned long length;
    char *sequence;
	NSString *name;
    NSString *description;
}

-(int) setSequence: (char*) init_seq;
-(int) setName: (NSString*) init_name;
-(int) setDescription: (NSString*) init_desc;

-(char*) getSequence;
-(NSString*) getName;
-(NSString*) getDescription;
-(unsigned long) getLength;

@end

@implementation DNASequence

-(int) setSequence:(char*)init_seq
{
	unsigned long len = strlen(init_seq);
	char b;
	sequence = calloc(sizeof(char), len+1);
	strncpy(sequence, init_seq, len);
    
	int i = 0;
	while (sequence[i] != '\0')
	{
		b = sequence[i];
		if (b != 'A' && b != 'T' && b != 'G' && b != 'C')
		{
			sequence = NULL;
			return -1;
		}
		i++;
	}

	length = len;
	return 0;
}


-(int) setName:(NSString *)init_name
{
    name = init_name;
    return 0;
    
}

-(int) setDescription:(NSString *)init_desc
{
    description = init_desc;
    return 0;
}

-(char*) getSequence
{
    return sequence;
}

-(NSString*) getName
{
    return name;
}

-(NSString*) getDescription
{
    return description;
}

-(unsigned long) getLength
{
    return length;
}

@end


@interface DuplexThermo: NSObject
{
    DNASequence *oligo1_seq;
    DNASequence *oligo2_seq;
    float temperature;
}
-(float) getDeltaG;
-(int) getThermoData: (float *) deltaS  deltaH:(float*) deltaH;


-(float) getTemperature;
-(int) setTemperature: (float) init_temp;
-(int) setOligo1: (DNASequence*) set_oligo;
-(int) setOligo2: (DNASequence*) set_oligo;
-(int) setOligo1withCString: (char*) set_oligo;
-(int) setOligo2withCString: (char*) set_oligo;


@end

@implementation DuplexThermo

-(float) getDeltaG
{
    //∆G = ∆H - T∆S
    float *deltaS = calloc(sizeof(float),1);
    float *deltaH = calloc(sizeof(float),1);
    [self getThermoData:deltaS deltaH:deltaH];

    
    return (*deltaH - ([self getTemperature]* (*deltaS) ));
    
    
}

// Finish adding non watson-crick pairs
-(int) getThermoData:(float *)deltaS deltaH:(float *)deltaH
{
    char *oligo1 = [oligo1_seq getSequence];
    char *oligo2 = [oligo2_seq getSequence];
    
    *deltaS = 0;
    *deltaH = 0;
    
    /*NSLog(@"oligo1:%s", oligo1);
    NSLog(@"oligo2:%s", oligo2);*/
    
    if (((strncmp(oligo1, "AC", 2) == 0)) || ((strncmp(oligo1, "AG", 2) == 0)) || ((strncmp(oligo1, "TG", 2) == 0)) || ((strncmp(oligo1, "TC", 2) == 0)) || ((strncmp(oligo1, "CC", 2) == 0)))
    {
        // Assuming that this is physically valid
        oligo1 = strrev(oligo1);
        oligo2 = strrev(oligo2);
    }
    
    if (((strncmp(oligo1, "TT", 2) == 0) && (strncmp(oligo2, "AA", 2)==0)) || ((strncmp(oligo1, "CC", 2) == 0) && (strncmp(oligo2, "GG", 2)==0)))
    {
        char *tmp = calloc(sizeof(char),3);
        tmp = oligo1;
        oligo1 = oligo2;
        oligo2 = tmp;
    }
    
    /*NSLog(@"oligo1 after swap/rev:%s", oligo1);
     NSLog(@"oligo2 after swa/rev:%s", oligo2);*/
    
    
    
    // Watson-Crick pairs
    /* From: Santa Lucia, J. (1998). A unified view of polymer, dumbbell, and
     oligonucleotide DNA nearest-neighbor thermodynamics . Proceedings of the
     National Academy of Science. */
    
    if ((strncmp(oligo1, "AA", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS =  -92.9;
        *deltaH = -33.1;
    }
    
    else if ((strncmp(oligo1, "AT", 2) == 0) && (strncmp(oligo2, "TA", 2) == 0))
    {
        *deltaS =  -85.4;
        *deltaH = -30.1;
    }
    
    else if ((strncmp(oligo1, "TA", 2) == 0) && (strncmp(oligo2, "AT", 2) == 0))
    {
        *deltaS =  -89.1;
        *deltaH = -30.1;
    }
    
    else if ((strncmp(oligo1, "CA", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS =  -95.0;
        *deltaH = -35.6;
    }
    
    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "CA", 2) == 0))
    {
        *deltaS =  -93.7;
        *deltaH = -35.1;
    }
    
    else if ((strncmp(oligo1, "CT", 2) == 0) && (strncmp(oligo2, "GA", 2) == 0))
    {
        *deltaS =  -87.9;
        *deltaH = -32.6;
    }
    
    else if ((strncmp(oligo1, "GA", 2) == 0) && (strncmp(oligo2, "CT", 2) == 0))
    {
        *deltaS =  -92.9;
        *deltaH = -34.3;
        
    }
    
    else if ((strncmp(oligo1, "CG", 2) == 0) && (strncmp(oligo2, "GC", 2) == 0))
    {
        *deltaS =  -113.8;
        *deltaH = -44.4;
    }
    
    else if ((strncmp(oligo1, "GC", 2) == 0) && (strncmp(oligo2, "CG", 2) == 0))
    {
        *deltaS =  -102.1;
        *deltaH = -41.0;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "CC", 2) == 0))
    {
        *deltaS =  -83.3;
        *deltaH = -33.5;
    }
    
    
    // Non Watson-Crick pairs
    /* From: Allawi, H. T., & SantaLucia, J. (1998). Thermodynamics of internal
     C.T mismatches in DNA. Nucleic Acids Research, 26(11), 2694–701.
     Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/10090733 */
    
    else if ((strncmp(oligo1, "AG", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS = 3.8;
        *deltaH = 4.2;
    }
    
    else if ((strncmp(oligo1, "AT", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = -35;
        *deltaH = -11;
    }
    
    else if ((strncmp(oligo1, "CG", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS = -49.0;
        *deltaH = -17;
        
    }
    
    else if ((strncmp(oligo1, "CT", 2) == 0) && (strncmp(oligo2, "GG", 2) == 0))
    {
        *deltaS = -34;
        *deltaH = -12;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "CT", 2) == 0))
    {
        *deltaS = 43.5;
        *deltaH = 14;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS = 68.2;
        *deltaH = 24;
    }
    
    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "CG", 2) == 0))
    {
        *deltaS = -51.5;
        *deltaH = -18;
    }
    
    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = 40;
        *deltaH = 17;
    }
    
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "AT", 2) == 0))
    {
        *deltaS = -7.1;
        *deltaH = -0.42;
    }
    
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS = -26;
        *deltaH = -5.9;
    }
    
    else if ((strncmp(oligo1, "TT", 2) == 0) && (strncmp(oligo2, "AG", 2) == 0))
    {
        *deltaS = -22;
        *deltaH = -5.4;
    }
    
    // More non-Watson-Crick sets
    /* From: Allawi, H. T., & SantaLucia, J. (1998). Nearest neighbor
     thermodynamic parameters for internal G.A mismatches in DNA.
     Biochemistry, 37(8), 2170–9. doi:10.1021/bi9724873 */
    
    else if ((strncmp(oligo1, "AA", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = -9.6;
        *deltaH = -2.5;
    }
    
    else if ((strncmp(oligo1, "AG", 2) == 0) && (strncmp(oligo2, "TA", 2) == 0))
    {
        *deltaS = -9.6;
        *deltaH = -2.9;
    }
    
    else if ((strncmp(oligo1, "CA", 2) == 0) && (strncmp(oligo2, "GG", 2) == 0))
    {
        *deltaS = -9.6;
        *deltaH = -2.9;
    }
    
    else if ((strncmp(oligo1, "CG", 2) == 0) && (strncmp(oligo2, "GA", 2) == 0))
    {
        *deltaS = -55.2;
        *deltaH = -16.7;
    }
    
    else if ((strncmp(oligo1, "GA", 2) == 0) && (strncmp(oligo2, "CG", 2) == 0))
    {
        *deltaS = -4.2;
        *deltaH = -2.5;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "CA", 2) == 0))
    {
        *deltaS = 13.4;
        *deltaH = 2.1;
    }
    
    else if ((strncmp(oligo1, "TA", 2) == 0) && (strncmp(oligo2, "AG", 2) == 0))
    {
        *deltaS = 2.9;
        *deltaH = 2.9;
    }
    
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "AA", 2) == 0))
    {
        *deltaS = 31;
        *deltaH = 13;
    }
    
    // More mismatches
    /* From: Allawi, H. T., & SantaLucia, J. (1997). Thermodynamics of internal
     G.T mismatches in DNA. Nucleic Acids Research, 36(11), 10581–10594.
     Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/10090733 */
    
    else if ((strncmp(oligo1, "AG", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS = 3.8;
        *deltaH = 4.2;
    }
    
    else if ((strncmp(oligo1, "AT", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = -35;
        *deltaH = -10;
    }
    
    else if ((strncmp(oligo1, "CG", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS = -49.0;
        *deltaH = -17.2;
    }
    
    else if ((strncmp(oligo1, "CT", 2) == 0) && (strncmp(oligo2, "GG", 2) == 0))
    {
        *deltaS = -33;
        *deltaH = -12;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "CT", 2) == 0))
    {
        *deltaS = 43.5;
        *deltaH = 14;
    }
    
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS = 68.2;
        *deltaH = 24;
    }
    
    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "CG", 2) == 0))
    {
        *deltaS = -51.5;
        *deltaH = -18;
    }
    
    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = 40;
        *deltaH = 17;
    }
    
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "AT", 2) == 0))
    {
        *deltaS = -7.1;
        *deltaH = -0.4;
    }
    
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS = -26;
        *deltaH = -5.9;
    }
    
    else if ((strncmp(oligo1, "TT", 2) == 0) && (strncmp(oligo2, "AG", 2) == 0))
    {
        *deltaS = -22;
        *deltaH = -5.4;
    }
    
    // More mismatches
    /* From: Allawi, H. T., & SantaLucia, J. (1998). Nearest-neighbor
     thermodynamics of internal A.C mismatches in DNA: sequence dependence
     and pH effects. Biochemistry, 37(26), 9435–44. doi:10.1021/bi9803729 */
    
    // These are for pH=7.0
    // TODO: Add in pH dependence? What about the others?
    
    else if ((strncmp(oligo1, "AA", 2) == 0) && (strncmp(oligo2, "TC", 2) == 0))
    {
        *deltaS = 19;
        *deltaH = 9.6;
    }
    
    else if ((strncmp(oligo1, "AC", 2) == 0) && (strncmp(oligo2, "TA", 2) == 0))
    {
        *deltaS = 61.1;
        *deltaH = 22;
    }
    
    else if ((strncmp(oligo1, "CA", 2) == 0) && (strncmp(oligo2, "GC", 2) == 0))
    {
        *deltaS = 15;
        *deltaH = 7.9;
    }
    
    else if ((strncmp(oligo1, "CC", 2) == 0) && (strncmp(oligo2, "GA", 2) == 0))
    {
        *deltaS = -2.5;
        *deltaH = 2.5;
    }
    
    else if ((strncmp(oligo1, "GA", 2) == 0) && (strncmp(oligo2, "CC", 2) == 0))
    {
        *deltaS = 59.4;
        *deltaH = 22;
    }
    
    else if ((strncmp(oligo1, "GC", 2) == 0) && (strncmp(oligo2, "CA", 2) == 0))
    {
        *deltaS = -16;
        *deltaH = -3;
    }
    
    else if ((strncmp(oligo1, "TA", 2) == 0) && (strncmp(oligo2, "AC", 2) == 0))
    {
        *deltaS = 33;
        *deltaH = 14;
    }
    
    else if ((strncmp(oligo1, "TC", 2) == 0) && (strncmp(oligo2, "AA", 2) == 0))
    {
        *deltaS = 84.5;
        *deltaH = 32;
    }
    
    
    /* More mismatches from: Peyret, N., Seneviratne, P. A., Allawi, H. T., & 
       Santalucia, J. (1999). Articles Nearest-Neighbor Thermodynamics and NMR 
       of DNA Sequences with Internal A.A, C.C, G.G, and T.T Mismatches. 
       Biochemistry, 38(12), 3468–3477. */
    
    // Assuming Table 2 is mislabeled and "∆H° (eu)" means "∆S° (eu)"
    
    else if ((strncmp(oligo1, "AA", 2) == 0) && (strncmp(oligo2, "TA", 2) == 0))
    {
        *deltaS = 7.1;
        *deltaH = 5.0;
    }
    
    else if ((strncmp(oligo1, "CA", 2) == 0) && (strncmp(oligo2, "GA", 2) == 0))
    {
        *deltaS = -18;
        *deltaH = -4;
    }

    else if ((strncmp(oligo1, "GA", 2) == 0) && (strncmp(oligo2, "CA", 2) == 0))
    {
        *deltaS = -41;
        *deltaH = -12;
    }
    
    else if ((strncmp(oligo1, "TA", 2) == 0) && (strncmp(oligo2, "AA", 2) == 0))
    {
        *deltaS = 54.0;
        *deltaH = 20;
    }

    else if ((strncmp(oligo1, "AC", 2) == 0) && (strncmp(oligo2, "TC", 2) == 0))
    {
        *deltaS = -18;
        *deltaH = 0.0001; // Error checks for zero, was pm 2.1
    }

    else if ((strncmp(oligo1, "CC", 2) == 0) && (strncmp(oligo2, "GC", 2) == 0))
    {
        *deltaS = -30;
        *deltaH = -6.3;
    }

    else if ((strncmp(oligo1, "GC", 2) == 0) && (strncmp(oligo2, "CC", 2) == 0))
    {
        *deltaS = 37;
        *deltaH = 15;
    }
             
    else if ((strncmp(oligo1, "TC", 2) == 0) && (strncmp(oligo2, "AC", 2) == 0))
    {
        *deltaS = 68.6;
        *deltaH = 26;
    }
             
    else if ((strncmp(oligo1, "AG", 2) == 0) && (strncmp(oligo2, "TG", 2) == 0))
    {
        *deltaS = -40;
        *deltaH = -13;
    }
             
    else if ((strncmp(oligo1, "CG", 2) == 0) && (strncmp(oligo2, "GG", 2) == 0))
    {
        *deltaS = -64.0;
        *deltaH = -21;
    }
             
    else if ((strncmp(oligo1, "GG", 2) == 0) && (strncmp(oligo2, "CG", 2) == 0))
    {
        *deltaS = -66.1;
        *deltaH = -25;
    }
             
    else if ((strncmp(oligo1, "TG", 2) == 0) && (strncmp(oligo2, "AG", 2) == 0))
    {
        *deltaS = 15;
        *deltaH = 6.7;
    }
             
    else if ((strncmp(oligo1, "AT", 2) == 0) && (strncmp(oligo2, "TT", 2) == 0))
    {
        *deltaS = -45.2;
        *deltaH = -11;
    }
             
    else if ((strncmp(oligo1, "CT", 2) == 0) && (strncmp(oligo2, "GT", 2) == 0))
    {
        *deltaS = -66.1;
        *deltaH = -21;
    }

    else if ((strncmp(oligo1, "GT", 2) == 0) && (strncmp(oligo2, "CT", 2) == 0))
    {
        *deltaS = -35;
        *deltaH = -9.2;
    }
             
    else if ((strncmp(oligo1, "TT", 2) == 0) && (strncmp(oligo2, "AT", 2) == 0))
    {
        *deltaS = -6.3;
        *deltaH = 0.8;
    }

    


        
    if ((*deltaS == 0) || (*deltaH==0))
    {
        //NSLog(@"Error: getThermoData:Duplex pair not found");
        return -1;
    }
    
    // Units above are in J/L*mol for ease of entry--convert to kJ/mol
    
    /*NSLog(@"deltaS:%f", *deltaS);
     NSLog(@"deltaH:%f", *deltaH);*/
    
    *deltaS = *deltaS/1000;
    
    return 0;
}

-(float) getTemperature
{
    return temperature;
}

-(int) setTemperature:(float)init_temp
{
    if (init_temp >= 0)
    {
        temperature = init_temp;
        return 0;
    }
    
    NSLog(@"setTemperature: Cannot set negative temperature--units are Kelvin");
    return -1;
}
-(int) setOligo1:(DNASequence *)set_oligo
{
    oligo1_seq = set_oligo;
    return 0;
}

-(int) setOligo2:(DNASequence *)set_oligo
{
    oligo2_seq = set_oligo;
    return 0;
}

-(int) setOligo1withCString:(char *)set_oligo
{
    DNASequence *tmpseq = [[DNASequence alloc] init];
    
    [tmpseq setSequence:set_oligo];
    oligo1_seq = tmpseq;
    
    return 0;
}

-(int) setOligo2withCString:(char *)set_oligo
{
    DNASequence *tmpseq = [[DNASequence alloc] init];
    
    [tmpseq setSequence:set_oligo];
    oligo2_seq = tmpseq;
    
    return 0;
}
@end


/*
 TODO:
 * Add support for dangling ends
 * Add support for long bulges
 * Add support for hairpins
 * Find the missing duplexes (AA/AT for example)
 * Fix end handling
 */

@interface OligoThermo: NSObject
{
    DNASequence *template;
    DNASequence *oligo;
    
    float temperature;
    
}
-(float) getDeltaGForPiece: (const char*) template_string;
-(float) getTemperature;

-(int) setTemperature: (float) init_temp;
-(int) setTemplate: (DNASequence*) init_template;
-(int) setOligo: (DNASequence*) init_oligo;

@end

@implementation OligoThermo


-(float) getLowestDeltaG
{
	float currentDeltaG = 0;
    float highestDeltaG = 0;
    
    
	char *template_string = calloc(sizeof(char), 10);
	char *oligo_string = calloc(sizeof(char), 10);
    
    template_string = [template getSequence];
    oligo_string = [oligo getSequence];
    
    unsigned long template_n = strlen(template_string);
    unsigned long oligo_n = strlen(oligo_string);
    
    char *binding_template = calloc(sizeof(char), oligo_n+1);
    
    for (int i = 0; i <= template_n - oligo_n; i++)
    {
        strncpy(binding_template, template_string+i, oligo_n);
		
        currentDeltaG = [self getDeltaGForPiece:binding_template];
        if (currentDeltaG < highestDeltaG)
            highestDeltaG = currentDeltaG;
        
    }
    return highestDeltaG;
}


-(float) getDeltaGForPiece: (const char*) template_string
{
    
    float deltaGtotal = 0;
    

    const char *oligo_string = [oligo getSequence];
    unsigned long oligo_n = strlen(oligo_string);
    
    DuplexThermo *duplexThermo = [[DuplexThermo alloc] init];
    [duplexThermo setTemperature:temperature];
    
    

    for (int j=0; j < oligo_n-1; j++)
    {
        char *seq1 = calloc(sizeof(char), 3);
        char *seq2 = calloc(sizeof(char), 3);
        
        strncpy(seq1, template_string+j, 2);
        strncpy(seq2, oligo_string+j, 2);
        
        [duplexThermo setOligo1withCString: seq1];
        [duplexThermo setOligo2withCString: seq2];
        
        float deltaGsingle = [duplexThermo getDeltaG];
        deltaGtotal = deltaGsingle + deltaGtotal;
        
        
        // account for initiation (first and last base)
        if ((j == 0) || (j == oligo_n-2))
        {
            if ((strncmp(seq1, "A", 1) == 0) || (strncmp(seq1, "T", 1)==0))
            {
                if ((strncmp(seq2, "A", 1) == 0) || (strncmp(seq2, "T", 1)==0))
                {
                    float deltaH = 9.62;
                    float deltaS = 17.15/1000;
                    float deltaGinit = deltaH - (temperature*deltaS);
                    deltaGtotal = deltaGtotal + deltaGinit;
                    
                    
                }
            }
            
            else if ((strncmp(seq1, "G", 1) == 0) || (strncmp(seq1, "C", 1)==0))
            {
                if ((strncmp(seq2, "G", 1) == 0) || (strncmp(seq2, "C", 1)==0))
                {
                    float deltaH = 0.4184;
                    float deltaS = -11.72/1000;
                    float deltaGinit = deltaH - (temperature*deltaS);
                    deltaGtotal = deltaGtotal + deltaGinit;
                }
            }
        }
        
        
    }
    

    
    return deltaGtotal;
}

-(float) getTemperature
{
    return temperature;
}

-(int) setTemperature:(float)init_temp
{
    if (init_temp >= 0)
    {
        temperature = init_temp;
        return 0;
    }
    
    NSLog(@"setTemperature: Error temperature cannot be less than one (units are K)");
    return -1;
}

-(int) setTemplate:(DNASequence *)init_template
{
    template = init_template;
    return 0;
}

-(int) setOligo:(DNASequence *)init_oligo
{
    oligo = init_oligo;
    return 0;
}

@end



char *open_sequence(int ID)
{
	char *seq = NULL;
	int len=0;
	
    
	// ID 0: pBSV2
	// ID 1: K12
    
	if (ID==0)
	{
		len = 7000;
	}
	else if (ID == 1)
	{
		len = 5000000;
	}
    
	seq = calloc(sizeof(char), len);
    
	FILE *fp;
	
	if (ID==0)
		fp = fopen("pBSV2", "r");
	else if (ID==1)
		fp = fopen("K12", "r");
    
	if (fp == NULL)
	{
		NSLog(@"Error opening sequence file");
		return NULL;
	}
    
	fgets(seq, len, fp);
    
	fclose(fp);
    
	return seq;
}


int main(int argc, const char * argv[])
{
	NSAutoreleasePool *pool = [NSAutoreleasePool new];
    
	printf("Starting program...\n");
	
	
	char *s1 = "CGTTCA";
	char *s2 = "GCAAGT";
	
    printf("Template sequence is:%s\n", s1);
    printf("Oligo sequence is:%s\n", s2);
    
	
    DNASequence *seq1 = [[DNASequence alloc] init];
    [seq1 setSequence:s1];
    
    DNASequence *seq2 = [[DNASequence alloc] init];
    [seq2 setSequence:s2];
    
    
    
    OligoThermo *oligthermo = [[OligoThermo alloc] init];
    [oligthermo setTemperature:310.15];
    [oligthermo setTemplate:seq1];
    [oligthermo setOligo:seq2];
    
    float deltaG = [oligthermo getLowestDeltaG];
    printf("Total deltaG:%f\n", deltaG);
    
    
    
    
    [pool drain];
    return 0;
}


