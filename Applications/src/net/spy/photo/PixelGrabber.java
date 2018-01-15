package net.spy.photo;

import java.io.File;

import de.amidar.AmidarSystem;
import de.amidar.cacheBench.CacheBenchParameters;

public class PixelGrabber {
	
	int[] valShort = {//five 8x8 blocks
		    -1,       -1,       -1,       -1,       -1,       -3090207, -4668463, -2959134, 
		    -2236182, -2236439, -2301976, -2367769, -2368025, -2433561, -2630683, -3090976,
		    -3288098, -2565147, -2499610, -3025439, -4143402, -3946280, -3420451, -3486244, 
		    -3617573, -3749158, -4406573, -4011817, -2565402, -2499610, -2762525, -3485476,
		    -3485732, -3354403, -3354403, -3354403, -3354403, -3354403, -3354403, -3354403, 
		    -4997683, -5391413, -920330,  -1,       -1,       -1,       -1,       -1,
		    -1,       -1,       -1,       -13740424,-5653816, -5390902, -14200717,-262915,
		    -1,       -1,       -1,       -1,       -1,       -8941145, -15778460,-14134668};/*,
		    -4733487, -4733487, -4733487, -4733487, -3418658, -1,       -1,       -11636595,
		    -15778460,-15778460,-15383960,-15252631,-15778460,-15778460,-15778460,-15712667,
		    -14989460,-15778460,-15778460,-9664352, -1,       -1,       -11636595,-15778460, 
		    -15778460,-15252631,-15186582,-15778460,-15778460,-6639938, -1,       -1,
		    -6179902, -15778460,-15712667,-12425595,-14266509,-15778460,-13806217,-723208,       
		    -1,       -1,       -11176559,-15778460,-15778460,-14726546,-15712667,-15778460,
		    -15778460,-7166024, -1,       -11176559,-15778460,-15778460,-14726546,-15712667, 
		    -15778460,-15778460,-7166024, -1,       -1,       -6179902, -15778460,-15712667,
		    -12425595,-14266509,-15778460,-13806217,-723208,  -1,       -11176559,-15778460, 
		    -15778460,-14726546,-15712667,-15778460,-15778460,-7166024, -1,       -1,
		    -1,       -1,       -1,       -1,       -1,       -1,       -1,       -1,       
		    -1,       -1,       -1,       -1,       -1,       -1,       -1,       -1,
		    -1,       -1,       -1,       -1,       -1,       -6706243, -2368024, -2104854, 
		    -2170646, -2236439, -2301976, -2367769, -2499611, -3025183, -3156513, -2499610,
		    -2499610, -3025439, -4208938, -3946280, -3420195, -3420451, -3420451, -3486244, 
		    -3617573, -3617829, -3617829, -3749158, -4340780, -4143146, -2565403, -2499610,
		    -2762524, -3485476, -3551268, -3354403, -3354403, -3354403, -3354403, -3354403, 
		    -3354403, -4734768, -5916987, -1,       -1,       -1,       -1,       -1,
		    -1,       -1,       -1,       -13477510,-6903109, -6705731, -14134668,-197123,       
		    -1,       -1,       -1,       -1,       -1,       -8941145, -15778460,-15778460,
		    -15778460,-15778460,-15778460,-15778460,-11373681,-1,       -1,       -11636595, 
		    -15778460,-12951424,-920330,  -854794,  -12688510,-15778460,-15252375,-2629659,
		    -131586,  -9466974, -15778460,-14858131,-1,       -1,       -11636595,-15778460, 
		    -13543302,-986123,  -591879,  -12491388,-15778460,-10913388,-1,       -65794,
		    -13937546,-15778460,-5982780, -1,       -657416,  -13411717,-15778460,-6508609,       
		    -1,       -197123,  -15449753,-15778460,-8283730, -65794,   -3352866, -15778460,
		    -15778460,-7166024, -197123,  -15449753,-15778460,-8283730, -65794,   -3352866, 
		    -15778460,-15778460,-7166024, -1,       -65794,   -13937546,-15778460,-5982780,
		    -1,       -657416,  -13411717,-15778460,-6508609, -197123,  -15449753,-15778460, 
		    -8283730, -65794,   -3352866, -15778460,-15778460,-7166024, -1,       -1,
		    -1,       -1,       -1,       -1,       -1,       -1,       -1,       -1,       
		    -1,       -1,       -1,       -1,       -1,       -1,       -1,       -1,	
		  
	}; //*/
	
	int xShort = 8;
	int yShort = 8;
	
	int scale = CacheBenchParameters.getBenchmarkScale();
	int factor = CacheBenchParameters.getBenchmarkScaleFactor();
	
	
	int [] valLong = AmidarSystem.readBMP(scale, factor);
	
	
	int xLong = AmidarSystem.getBMPwidth();
	int yLong = AmidarSystem.getBMPheight();
	
	int [] values;
	
	boolean longPic;

	public PixelGrabber(boolean longPic) {
		this.longPic = longPic;
	}

	public void grabPixels(int[] values){
		this.values = values;
		
		if(longPic){
			for(int i = 0; i<values.length; i++){
				this.values[i] = valLong[i];
			}
		} else {
			for(int i = 0; i<values.length; i++){
				this.values[i] = valShort[i];
			}
		}
	}
	
	public int getWidth(){
		if(longPic){
			return xLong;
		} else {
			return xShort;
		}
	}
	
	public int getHeight(){
		if(longPic){
			return yLong;
		} else {
			return yShort;
		}
	}
}
