package sa;

import cgramodel.CgraModel;

public interface CostFunction {

	
	public double getCost(CgraModel model);
	
	public boolean isValid();
}
