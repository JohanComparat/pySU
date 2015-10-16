
	def fit_2lines(self,wl,spec1d,err1d,a0=n.array([2795.,2802.]),lineName="AA",lineName2="BB", aTR=2797.5,cmin1=2650.,cmax1=2770.,cmin2=2807.,cmax2=2840.,prior=n.array([1.,1.,1e-18,-1e-18]),model="gaussian"):
		"""
		fits jointly two gaussian line profiles

		:param wl: wavelength (array, Angstrom)
		:param spec1d: flux observed in the broad band (array, f lambda)
		:param err1d: error on the flux observed in the broad band (array, f lambda)
		:param a0: expected position of the peak of the line in the observed frame (redshifted). 2 positions given.
		:param lineName: suffix characterizing the first line in the headers of the output
		:param lineName2: suffix characterizing the second line in the headers of the output
		:param DLC: wavelength extent to fit the continuum around the line. (def: 230 Angstrom)
		:param p0_sigma: prior on the line width in A (def: 15 A)
		:param p0_flux: prior on the line flux in erg/cm2/s/A (def: 8e-17)
		:param p0_share: prior on the share between the two [OII] lines. (def: 0.58)
		:param continuumSide: "left" = bluewards of the line or "right" = redwards of the line
		:param model: line model to be fitted : "gaussian", "lorentz" or "pseudoVoigt"

		Returns :
		 * array 1 with the parameters of the model
		 * array 2 with the model (wavelength, flux model)
		 * header corresponding to the array 1
		"""		
		header=" a0_"+lineName+" flux_"+lineName+" fluxErr_"+lineName+" sigma_"+lineName+" sigmaErr_"+lineName+" a0_"+lineName2+" flux_"+lineName2+" fluxErr_"+lineName2+" sigma_"+lineName2+" sigmaErr_"+lineName2+" continu_"+lineName+" continu_"+lineName2+" chi2_"+lineName+" ndof_"+lineName 

		outPutNF=n.array([a0[0], self.dV,self.dV,self.dV,self.dV, a0[1],self.dV,self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])

		domainLine=(wl>a0.min()-self.fitWidth)&(wl<a0.max()+self.fitWidth)
		domainContL=(wl>cmin1)&(wl<cmax1)
		domainContR=(wl>cmin2)&(wl<cmax2)
		if cmin1+50>wl.min() and cmax2-50<wl.max() and len(domainContL.nonzero()[0])>2 and len(domainContR.nonzero()[0])>2  :
			continuL=n.median(spec1d[domainContL])
			continuR=n.median(spec1d[domainContR])
			# model : absorption first, emission second
			def flG(aai,sig1,sig2,F1,F2): 
				aa=aai[(aai<aTR)]
				left=continuL + F1*(n.e**(-(aa-a0[0])**2./(2.*sig1**2.)))/(sig1*(2.*n.pi)**0.5)  
				aa=aai[(aai>=aTR)]
				right=continuR + F2*(n.e**(-(aa-a0[1])**2./(2.*sig2**2.)))/(sig2*(2.*n.pi)**0.5)
				return n.hstack((left, right))

			out = curve_fit(flG, wl[domainLine], spec1d[domainLine], p0=prior,sigma=err1d[domainLine],maxfev=500000000)
			if out[1].__class__==n.ndarray : # if the fit worked
				model1=flG(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3])
				var=err1d[domainLine]
				chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
				ndof=len(var)
				sig1=out[0][0]
				sig1Err=out[1][0][0]**0.5
				sig2=out[0][1]
				sig2Err=out[1][1][1]**0.5
				flux1=out[0][2]
				flux1Err=out[1][2][2]**0.5
				flux2=out[0][3]
				flux2Err=out[1][3][3]**0.5

				outPut=n.array([a0[0],flux1,flux1Err,sig1,sig1Err,a0[1],flux2,flux2Err,sig2,sig2Err,continuL,continuR,chi2,ndof])
				mod=n.array([wl[domainLine],model1])
								
				return outPut,mod,header
			else :
				return outPutNF,modNF,header
		else :
			return outPutNF,modNF,header


	def fit_doublet_position(wl,spec1d,err1d,aTR=4994.,a0=n.array([O3_4960,O3_5007]), cmin1=4400,cmax1=4800,cmin2=5025,cmax2=5150, p0=n.array([3.,3.,1e-16,3e-16,O3_4960,O3_5007]),model="gaussian"):
		""" re fits the [OIII] doublet to determine the centroid of the lines
		input single line fits of the [OIII] lines
		returns the position of the gaussians"""
		header=" a0_O3aDB flux_O3aDB fluxErr_O3aDB sigma_O3aDB sigmaErr_O3aDB a0_O3aDB_fit a0_O3aDB_fitErr a0_O3bDB flux_O3bDB fluxErr_O3bDB sigma_O3bDB sigmaErr_O3bDB a0_O3bDB_fit a0_O3bDB_fitErr cont_O3aDB cont_O3bDB chi2_O3DB"
		outPutNF=n.array([a0[0], self.dV,self.dV,self.dV, self.dV,self.dV,self.dV, a0[1],self.dV,self.dV,self.dV, self.dV,self.dV,self.dV, self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		#where the continuum and the lines are considered
		domainLine=(wl>a0.min()-self.fitWidth)&(wl<a0.max()+self.fitWidth)
		domainContL=(wl>cmin1)&(wl<cmax1)
		domainContR=(wl>cmin2)&(wl<cmax2)
		#print wl,spec1d,aTR,a0
		#print cmin1,cmax1,cmin2,cmax2
		#print spec1d[domainContL],spec1d[domainContR]
		# condition for fit
		if cmin1+50>wl.min() and cmax2-50<wl.max() and len(domainContL.nonzero()[0])>2 and len(domainContR.nonzero()[0])>2  :
			continuL=n.median(spec1d[domainContL])
			continuR=n.median(spec1d[domainContR])
			# model : absorption first, emission second
			def flG(aai,sig1,sig2,F1,F2,aa1,aa2): 
				aa=aai[(aai<aTR)]
				left=continuL + F1*(n.e**(-(aa-aa1)**2./(2.*sig1**2.)))/(sig1*(2.*n.pi)**0.5)  
				aa=aai[(aai>=aTR)]
				right=continuR + F2*(n.e**(-(aa-aa2)**2./(2.*sig2**2.)))/(sig2*(2.*n.pi)**0.5)
				return n.hstack((left, right))

			out = curve_fit(flG, wl[domainLine], spec1d[domainLine], p0=p0,sigma=err1d[domainLine],maxfev=500000000)
			# if the fit worked
			if out[1].__class__==n.ndarray :
				model1=flG(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5])
				var=err1d[domainLine]
				chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)/len(var)
				sig1=out[0][0]
				sig1Err=out[1][0][0]**0.5
				sig2=out[0][1]
				sig2Err=out[1][1][1]**0.5
				flux1=out[0][2]
				flux1Err=out[1][2][2]**0.5
				flux2=out[0][3]
				flux2Err=out[1][3][3]**0.5
				aa1=out[0][4]
				aa1Err=out[1][4][4]**0.5
				aa2=out[0][5]
				aa2Err=out[1][5][5]**0.5
				print aa1, aa2
				outPut=n.array([a0[0],flux1,flux1Err,sig1,sig1Err,aa1,aa1Err,a0[1],flux2,flux2Err,sig2,sig2Err,aa2,aa2Err,continuL,continuR,chi2])
				mod=n.array([wl[domainLine],model1])
				return outPut,mod, header
			else :
				return outPutNF,modNF,header
		else :
			return outPutNF,modNF,header


	def fit_3lines(wl,spec1d,err1d,a0=n.array([6564.61,6549.859,6585.268]),lineName="AAA",lineName2="BBB",lineName3="CCC",aTR=6575.,cmin1=6400.,cmax1=6525.,cmin2=6600.,cmax2=6700.,prior=n.array([2.,4.,2.,1e-17,1e-16,1e-17]),model="gaussian"):
		header=" a0_"+lineName+" flux_"+lineName+" fluxErr_"+lineName+" sigma_"+lineName+" sigmaErr_"+lineName+" a0_"+lineName2+" flux_"+lineName2+" fluxErr_"+lineName2+" sigma_"+lineName2+" sigmaErr_"+lineName2+" a0_"+lineName3+" flux_"+lineName3+" fluxErr_"+lineName3+" sigma_"+lineName3+" sigmaErr_"+lineName3+" continu_"+lineName+" continu_"+lineName3+" chi2_"+lineName+" ndof_"+lineName
		outPutNF=n.array([a0[0],self.dV,self.dV,self.dV,self.dV, a0[1],self.dV,self.dV,self.dV,self.dV, a0[2],self.dV,self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])
		domainLine=(wl>a0.min()-self.fitWidth)&(wl<a0.max()+self.fitWidth)
		domainContL=(wl>cmin1)&(wl<cmax1)
		domainContR=(wl>cmin2)&(wl<cmax2)
		if cmin1+50>wl.min() and cmax2-50<wl.max() and len(domainContL.nonzero()[0])>2 and len(domainContR.nonzero()[0])>2  :
			continuL=n.median(spec1d[domainContL])
			continuR=n.median(spec1d[domainContR])
			# model : absorption first, emission second
			def flG(aai,sig1,sig2,sig3,F1,F2,F3 ): 
				aa=aai[(aai<aTR)]
				left=continuL + F1*(n.e**(-(aa-a0[0])**2./(2.*sig1**2.)))/(sig1*(2.*n.pi)**0.5) + F2*(n.e**(-(aa-a0[1])**2./(2.*sig2**2.)))/(sig2*(2.*n.pi)**0.5) 
				aa=aai[(aai>=aTR)]
				right=continuR + F3*(n.e**(-(aa-a0[2])**2./(2.*sig3**2.)))/(sig3*(2.*n.pi)**0.5) 
				return n.hstack((left, right))

			out = curve_fit(flG, wl[domainLine], spec1d[domainLine], p0=prior,sigma=err1d[domainLine],maxfev=500000000)
			if out[1].__class__==n.ndarray : # if the fit worked
				model1=flG(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5])
				var=err1d[domainLine]
				chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
							ndof=len(var)
				sig1=out[0][0]
				sig1Err=out[1][0][0]**0.5
				sig2=out[0][1]
				sig2Err=out[1][1][1]**0.5
				sig3=out[0][2]
				sig3Err=out[1][2][2]**0.5
				flux1=out[0][3]
				flux1Err=out[1][3][3]**0.5
				flux2=out[0][4]
				flux2Err=out[1][4][4]**0.5
				flux3=out[0][5]
				flux3Err=out[1][5][5]**0.5
				outPut=n.array([a0[0],flux1,flux1Err,sig1,sig1Err,a0[1],flux2,flux2Err,sig2,sig2Err,a0[2],flux3,flux3Err,sig3,sig3Err ,continuL,continuR,chi2,ndof])
				mod=n.array([wl[domainLine],model1])
				return outPut,mod,header
			else :
				return outPutNF,modNF,header
		else :
			return outPutNF,modNF,header


	def fit_4lines(wl,spec1d,err1d,a0=n.array([2586.1,2599.6,2612.5,2626.3]),lineName="AAAA",lineName2="BBBB",lineName3="CCCC",lineName4="DDDD",aTR=2606.,cmin1=2400.,cmax1=2550.,cmin2=2650.,cmax2=2770.,prior=n.array([1.,1.,1.,1.,1e-18,-1e-18,1e-18,-1e-18]),model="gaussian"):
		header=" a0_"+lineName+" flux_"+lineName+" fluxErr_"+lineName+" sigma_"+lineName+" sigmaErr_"+lineName+" a0_"+lineName2+" flux_"+lineName2+" fluxErr_"+lineName2+" sigma_"+lineName2+" sigmaErr_"+lineName2+" a0_"+lineName3+" flux_"+lineName3+" fluxErr_"+lineName3+" sigma_"+lineName3+" sigmaErr_"+lineName3+" a0_"+lineName4+" flux_"+lineName4+" fluxErr_"+lineName4+" sigma_"+lineName4+" sigmaErr_"+lineName4+" continu_"+lineName+" continu_"+lineName4+" chi2_"+lineName+" ndof_"+lineName

		outPutNF=n.array([a0[0],self.dV,self.dV,self.dV,self.dV, a0[1],self.dV,self.dV,self.dV,self.dV, a0[2],self.dV,self.dV,self.dV,self.dV, a0[3],self.dV,self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])

		domainLine=(wl>a0.min()-self.fitWidth)&(wl<a0.max()+self.fitWidth)
		domainContL=(wl>cmin1)&(wl<cmax1)
		domainContR=(wl>cmin2)&(wl<cmax2)
		if cmin1+50>wl.min() and cmax2-50<wl.max() and len(domainContL.nonzero()[0])>2 and len(domainContR.nonzero()[0])>2  :
			continuL=n.median(spec1d[domainContL])
			continuR=n.median(spec1d[domainContR])
			# model : absorption first, emission second
			def flG(aai,sig1,sig2,sig3,sig4,F1,F2,F3,F4 ): 
				aa=aai[(aai<aTR)]
				left=continuL + F1*(n.e**(-(aa-a0[0])**2./(2.*sig1**2.)))/(sig1*(2.*n.pi)**0.5) + F2*(n.e**(-(aa-a0[1])**2./(2.*sig2**2.)))/(sig2*(2.*n.pi)**0.5) 
				aa=aai[(aai>=aTR)]
				right=continuR + F3*(n.e**(-(aa-a0[2])**2./(2.*sig3**2.)))/(sig3*(2.*n.pi)**0.5) + F4*(n.e**(-(aa-a0[3])**2./(2.*sig4**2.)))/(sig4*(2.*n.pi)**0.5) 
				return n.hstack((left, right))

			out = curve_fit(flG, wl[domainLine], spec1d[domainLine], p0=prior,sigma=err1d[domainLine],maxfev=500000000)
			if out[1].__class__==n.ndarray : # if the fit worked
				model1=flG(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5],out[0][6],out[0][7])
				var=err1d[domainLine]
				chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
							ndof=len(var)
				sig1=out[0][0]
				sig1Err=out[1][0][0]**0.5
				sig2=out[0][1]
				sig2Err=out[1][1][1]**0.5
				sig3=out[0][2]
				sig3Err=out[1][2][2]**0.5
				sig4=out[0][3]
				sig4Err=out[1][3][3]**0.5
				flux1=out[0][4]
				flux1Err=out[1][4][4]**0.5
				flux2=out[0][5]
				flux2Err=out[1][5][5]**0.5
				flux3=out[0][6]
				flux3Err=out[1][6][6]**0.5
				flux4=out[0][7]
				flux4Err=out[1][7][7]**0.5
				outPut=n.array([a0[0],flux1,flux1Err,sig1,sig1Err,a0[1],flux2,flux2Err,sig2,sig2Err,a0[2],flux3,flux3Err,sig3,sig3Err,a0[3],flux4,flux4Err,sig4,sig4Err ,continuL,continuR,chi2,ndof])
				mod=n.array([wl[domainLine],model1])
				return outPut,mod,header
			else :
				return outPutNF,modNF,header
		else :
			return outPutNF,modNF,header

	def fit_6lines(wl,spec1d,err1d,a0=n.array([2326.7,2343.7,2365.3,2374.3,2382.2, 2396.2]),lineName="AAAA",lineName2="BBBB", lineName3="CCCC",lineName4="DDDD", lineName5="EEEE",lineName6="FFFF",aTR=2370.,cmin1=2080.,cmax1=2240.,cmin2=2400.,cmax2=2550.,prior=n.array([1.,1.,1.,1.,1.,1e-18,-1e-18,1e-18,-1e-18,-1e-18,1e-18,1.]),model="gaussian"):
		header=" a0_"+lineName+" flux_"+lineName+" fluxErr_"+lineName+" sigma_"+lineName+" sigmaErr_"+lineName+" a0_"+lineName2+" flux_"+lineName2+" fluxErr_"+lineName2+" sigma_"+lineName2+" sigmaErr_"+lineName2+" a0_"+lineName3+" flux_"+lineName3+" fluxErr_"+lineName3+" sigma_"+lineName3+" sigmaErr_"+lineName3+" a0_"+lineName4+" flux_"+lineName4+" fluxErr_"+lineName4+" sigma_"+lineName4+" sigmaErr_"+lineName4+" a0_"+lineName5+" flux_"+lineName5+" fluxErr_"+lineName5+" sigma_"+lineName5+" sigmaErr_"+lineName5+" a0_"+lineName6+" flux_"+lineName6+" fluxErr_"+lineName6+" sigma_"+lineName6+" sigmaErr_"+lineName6+" continu_"+lineName+" continu_"+lineName6+" chi2_"+lineName+" ndof_"+lineName
		outPutNF=n.array([a0[0],self.dV,self.dV,self.dV,self.dV, a0[1],self.dV,self.dV,self.dV,self.dV, a0[2],self.dV,self.dV,self.dV,self.dV, a0[3],self.dV,self.dV,self.dV,self.dV, a0[4],self.dV,self.dV,self.dV,self.dV, a0[5],self.dV,self.dV,self.dV,self.dV, self.dV,self.dV,self.dV,self.dV])
		modNF=n.array([self.dV,self.dV])

		domainLine=(wl>a0.min()-self.fitWidth)&(wl<a0.max()+self.fitWidth)
		domainContL=(wl>cmin1)&(wl<cmax1)
		domainContR=(wl>cmin2)&(wl<cmax2)
		if cmin1+50>wl.min() and cmax2-50<wl.max() and len(domainContL.nonzero()[0])>2 and len(domainContR.nonzero()[0])>2  :
			continuL=n.median(spec1d[domainContL])
			continuR=n.median(spec1d[domainContR])
			# model : absorption first, emission second
			def flG(aai,sig1,sig2,sig3,sig4,sig5,F1,F2,F3,F4,F5,F6,sig6 ): 
				aa=aai[(aai<aTR)]
				left=continuL + F1*(n.e**(-(aa-a0[0])**2./(2.*sig1**2.)))/(sig1*(2.*n.pi)**0.5) + F2*(n.e**(-(aa-a0[1])**2./(2.*sig2**2.)))/(sig2*(2.*n.pi)**0.5) + F3*(n.e**(-(aa-a0[2])**2./(2.*sig3**2.)))/(sig3*(2.*n.pi)**0.5)
				aa=aai[(aai>=aTR)]
				right=continuR + F4*(n.e**(-(aa-a0[3])**2./(2.*sig4**2.)))/(sig4*(2.*n.pi)**0.5) + F5*(n.e**(-(aa-a0[4])**2./(2.*sig5**2.)))/(sig5*(2.*n.pi)**0.5) + F6*(n.e**(-(aa-a0[5])**2./(2.*sig6*2.)))/(sig6*(2.*n.pi)**0.5)
				return n.hstack((left, right))

			out = curve_fit(flG, wl[domainLine], spec1d[domainLine], p0=prior,sigma=err1d[domainLine],maxfev=500000000)
			if out[1].__class__==n.ndarray : # if the fit worked
				model1=flG(wl[domainLine],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4],out[0][5],out[0][6],out[0][7],out[0][8],out[0][9],out[0][10],out[0][11])
				var=err1d[domainLine]
				chi2=n.sum(abs(model1-spec1d[domainLine])**2./var**2.)
				ndof=len(var)
				sig1=out[0][0]
				sig1Err=out[1][0][0]**0.5
				sig2=out[0][1]
				sig2Err=out[1][1][1]**0.5
				sig3=out[0][2]
				sig3Err=out[1][2][2]**0.5
				sig4=out[0][3]
				sig4Err=out[1][3][3]**0.5
				sig5=out[0][4]
				sig5Err=out[1][4][4]**0.5
				flux1=out[0][5]
				flux1Err=out[1][5][5]**0.5
				flux2=out[0][6]
				flux2Err=out[1][6][6]**0.5
				flux3=out[0][7]
				flux3Err=out[1][7][7]**0.5
				flux4=out[0][8]
				flux4Err=out[1][8][8]**0.5
				flux5=out[0][9]
				flux5Err=out[1][9][9]**0.5
				flux6=out[0][10]
				flux6Err=out[1][10][10]**0.5
				sig6=out[0][11]
				sig6Err=out[1][11][11]**0.5

				outPut=n.array([a0[0],flux1,flux1Err,sig1,sig1Err,a0[1],flux2,flux2Err,sig2,sig2Err,a0[2],flux3,flux3Err,sig3,sig3Err,a0[3],flux4,flux4Err,sig4,sig4Err,a0[4],flux5,flux5Err,sig5,sig5Err,a0[5],flux6,flux6Err,sig6,sig6Err ,continuL,continuR,chi2,ndof])
				mod=n.array([wl[domainLine],model1])

				return outPut,mod,header
			else :
				return outPutNF,modNF,header
		else :
			return outPutNF,modNF,header



	def fit_4000D(wl,spec1d,err1d,z,intLim4k=n.array([3600-150, 3600, 4140, 4140+150])):
		""" returns the log of the ratio of the integrated fluxes
		 and the log of the ratio of the flux median density"""
		header=" logd4k logd4kU logd4kL" #left4kMedian left4kMedianU left4kMedianL right4kMedian right4kMedianU right4kMedianL "
		modNF=n.array([self.dV,self.dV])
		outPutNF=n.array([self.dV,self.dV,self.dV])
		if wl.min()<intLim4k[0]*(1+z) and wl.max()>intLim4k[3]*(1+z):
			toInt=interp1d(wl,spec1d)
			left=quad(toInt,intLim4k[0]*(1+z),intLim4k[1]*(1+z))[0]
			right=quad(toInt,intLim4k[2]*(1+z),intLim4k[3]*(1+z))[0]

			toInt=interp1d(wl,spec1d+err1d)
			leftU=quad(toInt,intLim4k[0]*(1+z),intLim4k[1]*(1+z))[0]
			rightU=quad(toInt,intLim4k[2]*(1+z),intLim4k[3]*(1+z))[0]

			toInt=interp1d(wl,spec1d-err1d)
			leftL=quad(toInt,intLim4k[0]*(1+z),intLim4k[1]*(1+z))[0]
			rightL=quad(toInt,intLim4k[2]*(1+z),intLim4k[3]*(1+z))[0]

			d4=n.log10(right/left)
			d4U=n.log10(rightL/leftU)
			d4L=n.log10(rightU/leftL)

			leftMedian=n.median(spec1d[(wl>intLim4k[0]*(1+z))&(wl<intLim4k[1]*(1+z))])
			leftMedianU=n.median(spec1d[(wl>intLim4k[0]*(1+z))&(wl<intLim4k[1]*(1+z))]+err1d[(wl>intLim4k[0]*(1+z))&(wl<intLim4k[1]*(1+z))])
			leftMedianL=n.median(spec1d[(wl>intLim4k[0]*(1+z))&(wl<intLim4k[1]*(1+z))]-err1d[(wl>intLim4k[0]*(1+z))&(wl<intLim4k[1]*(1+z))])

			rightMedian=n.median(spec1d[(wl>intLim4k[2]*(1+z))&(wl<intLim4k[3]*(1+z))])
			rightMedianU=n.median(spec1d[(wl>intLim4k[2]*(1+z))&(wl<intLim4k[3]*(1+z))]+err1d[(wl>intLim4k[2]*(1+z))&(wl<intLim4k[3]*(1+z))])
			rightMedianL=n.median(spec1d[(wl>intLim4k[2]*(1+z))&(wl<intLim4k[3]*(1+z))]-err1d[(wl>intLim4k[2]*(1+z))&(wl<intLim4k[3]*(1+z))])

			dM4=n.log10(rightMedian/leftMedian)
			dM4U=n.log10(rightMedianL/leftMedianU)
			dM4L=n.log10(rightMedianU/leftMedianL)
	
			outPut=n.array([d4,d4U,d4L]) #,leftMedian,leftMedianU,leftMedianL,rightMedian,rightMedianU,rightMedianL
			mod=n.array([self.dV,self.dV])
			return outPut,mod,header
		else:
			return outPutNF,modNF,header


	def fit_UV_slope(wl,spec1d,err1d,z,intLimUV=n.array([2000,2200,3000,3200,3400,3600,4100, 4300,4500,4700])):
		""" returns the luminosities in UV wavelength bands"""
		header=" L2100 L2100Err L3100 L3100Err L3500 L3500Err L4200 L4200Err L4600 L4600Err"
		modNF=n.array([self.dV,self.dV])

			if wl.min()<intLimUV[0]*(1+z) and wl.max()>intLimUV[1]*(1+z) :
					toInt=interp1d(wl,spec1d)
					L2100=quad(toInt,intLimUV[0]*(1+z),intLimUV[1]*(1+z))[0]
					toInt=interp1d(wl,spec1d+err1d)
					L2100U=quad(toInt,intLimUV[0]*(1+z),intLimUV[1]*(1+z))[0]
					toInt=interp1d(wl,spec1d-err1d)
					L2100L=quad(toInt,intLimUV[0]*(1+z),intLimUV[1]*(1+z))[0]
			else:
					L2100,L2100U,L2100L=self.dV,self.dV,self.dV

			if wl.min()<intLimUV[2]*(1+z) and wl.max()>intLimUV[3]*(1+z) :
					toInt=interp1d(wl,spec1d)
					L3100=quad(toInt,intLimUV[2]*(1+z),intLimUV[3]*(1+z))[0]
					toInt=interp1d(wl,spec1d+err1d)
					L3100U=quad(toInt,intLimUV[2]*(1+z),intLimUV[3]*(1+z))[0]
					toInt=interp1d(wl,spec1d-err1d)
					L3100L=quad(toInt,intLimUV[2]*(1+z),intLimUV[3]*(1+z))[0]
			else:
					L3100,L3100U,L3100L=self.dV,self.dV,self.dV

			if wl.min()<intLimUV[4]*(1+z) and wl.max()>intLimUV[5]*(1+z) :
					toInt=interp1d(wl,spec1d)
					L3500=quad(toInt,intLimUV[4]*(1+z),intLimUV[5]*(1+z))[0]
					toInt=interp1d(wl,spec1d+err1d)
					L3500U=quad(toInt,intLimUV[4]*(1+z),intLimUV[5]*(1+z))[0]
					toInt=interp1d(wl,spec1d-err1d)
					L3500L=quad(toInt,intLimUV[4]*(1+z),intLimUV[5]*(1+z))[0]
			else:
					L3500,L3500U,L3500L=self.dV,self.dV,self.dV

			if wl.min()<intLimUV[6]*(1+z) and wl.max()>intLimUV[7]*(1+z) :
					toInt=interp1d(wl,spec1d)
					L4200=quad(toInt,intLimUV[6]*(1+z),intLimUV[7]*(1+z))[0]
					toInt=interp1d(wl,spec1d+err1d)
					L4200U=quad(toInt,intLimUV[6]*(1+z),intLimUV[7]*(1+z))[0]
					toInt=interp1d(wl,spec1d-err1d)
					L4200L=quad(toInt,intLimUV[6]*(1+z),intLimUV[7]*(1+z))[0]
			else:
					L4200,L4200U,L4200L=self.dV,self.dV,self.dV

			if wl.min()<intLimUV[8]*(1+z) and wl.max()>intLimUV[9]*(1+z) :
					toInt=interp1d(wl,spec1d)
					L4600=quad(toInt,intLimUV[8]*(1+z),intLimUV[9]*(1+z))[0]
					toInt=interp1d(wl,spec1d+err1d)
					L4600U=quad(toInt,intLimUV[8]*(1+z),intLimUV[9]*(1+z))[0]
					toInt=interp1d(wl,spec1d-err1d)
					L4600L=quad(toInt,intLimUV[8]*(1+z),intLimUV[9]*(1+z))[0]
			else:
					L4600,L4600U,L4600L=self.dV,self.dV,self.dV
	
		outPut=n.array([L2100, (L2100U-L2100L)/2., L3100, (L3100U-L3100L)/2., L3500, (L3500U-L3500L)/2., L4200, (L4200U-L4200L)/2., L4600, (L4600U-L4600L)/2.])
		return outPut,modNF,header	


