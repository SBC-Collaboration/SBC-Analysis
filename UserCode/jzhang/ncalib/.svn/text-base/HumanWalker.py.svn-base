import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QHBoxLayout, QSizePolicy, QMessageBox, \
    QWidget, QPushButton, QCheckBox, QTableWidget, QDoubleSpinBox, QSpinBox, QFileDialog
from PyQt5.QtCore import QSize, QRect, pyqtSignal
from PyQt5.QtGui import QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import numpy as np
import random
import pickle
# import ipdb
import SBCcode.MonteCarlo.PICOcalGlobalLikelihood_v2_jz as pcgl


class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.left = 10
        self.top = 30
        self.title = 'PICO Nuclear Calibration'
        self.width = 1200
        self.height = 800
        self.chain = []
        self.chi2 = []

        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setWindowIcon(QIcon("folder.png"))

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        # self.m = PlotCanvas(self, width=5, height=4)
        # self.m.move(0,0)
        # self.m.plot([1,3],[[3,4],[5,6]],'hello')

        # self.m1 = PlotCanvas(self, width=5, height=4)
        # self.m1.move(0,400)

        self.figure = Figure()
        self.ax_eff = self.figure.add_subplot(211)
        self.ax_chi2 = self.figure.add_subplot(212)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self.central_widget)

        left_box = QVBoxLayout()
        left_box.addWidget(self.toolbar)
        left_box.addWidget(self.canvas)

        table_holder = QWidget()
        self.table = ParamTable(table_holder)

        # # Just some button connected to `plot` method
        # self.button = QPushButton('Plot', self.central_widget)
        # self.button.clicked.connect(self.plot)

        # layout.addWidget(self.button)
        # layout.setGeometry(QRect(0, 0, 1200, 800))
        # layout.setGeometry(QRect(0, 0, 600, 800))

        param_box = QHBoxLayout()

        self.step_box = QSpinBox(self)
        self.step_box.setToolTip('Show parameters at step #')
        self.step_box.setPrefix('Step: ')
        self.step_box.setKeyboardTracking(False)
        self.step_box.valueChanged.connect(self.show_step)
        self.step_box.setMaximum(0)
        param_box.addWidget(self.step_box)

        button = QPushButton('Eff', self)
        button.setToolTip('Update Efficiency Plot')
        button.clicked.connect(self.plot_eff)
        param_box.addWidget(button)

        button = QPushButton('Chi^2', self)
        button.setToolTip('Update Chi^2 Plot')
        button.clicked.connect(self.update_plots)
        param_box.addWidget(button)

        self.is_auto = QCheckBox(self)
        self.is_auto.setText('Auto Update')
        self.is_auto.setToolTip('Automatically Updates Plots')
        self.is_auto.stateChanged.connect(self.auto_update)
        param_box.addWidget(self.is_auto)

        button_box = QHBoxLayout()

        button = QPushButton('Add Node', self)
        button.setToolTip('Add a (Er, P) node')

        button.clicked.connect(self.table.add_node)
        button_box.addWidget(button)

        button = QPushButton('Remove Node', self)
        button.setToolTip('Remove an Efficiency')
        button.clicked.connect(self.table.remove_node)
        button_box.addWidget(button)

        button = QPushButton('Minuit', self)
        button.setToolTip('Run Minuit')
        button.clicked.connect(self.run_minuit)
        # button.move(700,600)
        button_box.addWidget(button)

        button = QPushButton('Save', self)
        button.setToolTip('Saves the history')
        button.clicked.connect(self.save_history)
        button_box.addWidget(button)

        button = QPushButton('Load', self)
        button.setToolTip('Load saved history')
        button.clicked.connect(self.load_history)
        button_box.addWidget(button)

        right_box = QVBoxLayout()
        right_box.addWidget(table_holder)
        right_box.addLayout(param_box)
        right_box.addLayout(button_box)

        layout = QHBoxLayout()
        layout.addLayout(left_box)
        layout.addLayout(right_box)

        self.central_widget.setLayout(layout)

        self.show()

    # def on_changed_value(self, value):
    #     self.m.plot()
    #     self.table.on_changed_value()
    #     print(value)
    def show_step(self, i_step):
        self.table.set_value(self.chain[i_step])
        self.plot_eff()

    def auto_update(self, state):
        if self.is_auto.isChecked():
            self.table.value_changed.connect(self.update_plots)
        else:
            self.table.value_changed.disconnect(self.update_plots)

    def update_plots(self):
        self.get_chi2()
        self.plot_eff()
        self.plot_chi2()

    def plot_eff(self):
        n_Eff = self.table.n_Eff
        n_col = self.table.n_col
        values = np.array(self.table.get_value())
        Eff_ix = np.arange(n_Eff) * n_col
        Eff = values[Eff_ix]
        Er_C = values[Eff_ix + 1]
        Er_F = values[Eff_ix + 2]

        self.ax_eff.clear()
        self.ax_eff.plot(Er_F, Eff, 'r*-')
        self.ax_eff.plot(Er_C, Eff, 'b*-')
        self.ax_eff.legend(('F', 'C'))
        self.ax_eff.set_xlabel('Er (keV)')
        self.ax_eff.set_ylabel('Efficiency')
        self.ax_eff.set_title('Bubble Nucleation Efficiency')

        self.canvas.draw()

    def plot_chi2(self):

        self.ax_chi2.clear()
        self.ax_chi2.plot(self.chi2, 'r*-')
        self.ax_chi2.legend(('Chi^2',))
        self.ax_chi2.set_xlabel('Step #')
        self.ax_chi2.set_ylabel('Chi^2')
        # self.ax_chi2.set_title('Bubble Nucleation Efficiency')

        self.canvas.draw()

    def get_chi2(self):
        # convert gui displayed values to theta in the likelihood
        n_Eff = self.table.n_Eff
        n_col = self.table.n_col
        values = np.array(self.table.get_value())
        Eff_ix = np.arange(n_Eff) * n_col

        pcgl.p_fenceposts = values[Eff_ix]

        cut = np.ones(values.shape, dtype=np.bool)
        cut[Eff_ix] = False
        v1_theta = values[cut]

        chi2 = -2 * pcgl.PICOcalLL(self.v2_theta(v1_theta))

        if not (chi2 == np.inf or chi2 == -np.inf):
            self.chi2.append(chi2)
            self.chain.append(values)
            self.step_box.setMaximum(len(self.chi2) - 1)

        print(chi2)
        return chi2

    # convert v1 theta to v2 theta in PICO likelihood
    def v2_theta(self, v1_theta):
        var = v1_theta.copy()
        n_Epts = self.table.n_Eff * (self.table.n_col - 1)
        pcgl.n_Epts = n_Epts  # need to update this
        Epts = np.reshape(var[:n_Epts], (pcgl.p_fenceposts.size,
                                         pcgl.threshold_fenceposts.size,
                                         pcgl.species_list.size))

        dEpts = np.zeros(Epts.shape, dtype=Epts.dtype)
        dEpts[1:, :, :] = np.diff(Epts, 1, axis=0)
        dEpts[0, :, :] = Epts[0, :, :]

        var[:n_Epts] = np.reshape(dEpts, (pcgl.p_fenceposts.size *
                                          pcgl.threshold_fenceposts.size *
                                          pcgl.species_list.size))
        print(var)
        return var

    def save_history(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "QFileDialog.getSaveFileName()", "",
                                                  "All Files (*);;Pickle Files (*.p)", options=options)

        if fileName:
            print('Saving history to ' + fileName)
            with open(fileName, "wb") as pickling_on:
                pickle.dump((self.chain, self.chi2), pickling_on)
                pickling_on.close()
        # save list of walkers = chain, chi^2, likelihood, nu's

    def load_history(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            print('Loading history from ' + fileName)
            with open(fileName, "rb") as pickle_off:
                self.chain, self.chi2 = pickle.load(pickle_off)

        if len(self.chain) > 0:
            self.table.set_value(self.chain[-1])
            self.plot_eff()
            self.plot_chi2()

    def run_minuit(self):
        print('minuit')


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        self.toolbar = NavigationToolbar(self, parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # self.plot()

    def plot(self, x, y, title):
        data = [random.random() for i in range(25)]
        ax = self.figure.add_subplot(111)
        ax.plot(x, y, 'r-')
        ax.set_title(title)
        self.draw()


class ParamTable(QTableWidget):
    value_changed = pyqtSignal(float)

    def __init__(self, parent=None, n_Eff=5, n_nuisance=14,
                 values=[0.0, 8.1508, 3.8139,
                         0.2, 9.3993, 3.8495,
                         0.5, 12.5022, 5.1765,
                         0.8, 14.4248, 6.4186,
                         1.0, 22.0866, 12.5607,
                         -0.2518, -2.3705, 2.1569, -0.6839, -1.9952, -0.4400, 1.7738, -3.5711, 1.2013, -1.1032, 1.3775,
                         -1.8020, -0.1159, 2.3455]):

        super().__init__(parent)
        self.n_Eff = n_Eff
        self.n_nuisance = n_nuisance
        self.n_col = 3

        self.tableWidget = QTableWidget(parent)
        # self.tableWidget.setGeometry(QRect(0, 0, 100 * (self.n_col + 1), 30 * (self.n_row() + 1)))
        self.tableWidget.setColumnCount(self.n_col)
        # self.tableWidget.setRowCount(self.n_row())
        self.tableWidget.setRowCount(0)
        self.tableWidget.setObjectName("tableWidget")
        self.fill_table()
        self.set_value(values)
        self.connect_signal()
        self.resize_table()

    def n_row(self):
        return self.n_Eff + np.ceil(self.n_nuisance / self.n_col).astype('int') + 2

    def n_params(self):
        return self.n_Eff * self.n_col + self.n_nuisance

    def fill_table(self):
        for i_row in range(self.n_row()):
            self.add_row(i_row)

    def fill_table2(self, values=[]):

        if len(values) == 0:
            values = [0.0] * self.n_params()

        self.tableWidget.setColumnCount(self.n_col)
        self.tableWidget.setRowCount(self.n_row())

        # column header
        item = QTableWidgetItem()
        item.setText('Eff')
        self.tableWidget.setHorizontalHeaderItem(0, item)

        item = QTableWidgetItem()
        item.setText('Er_F')
        self.tableWidget.setHorizontalHeaderItem(1, item)

        item = QTableWidgetItem()
        item.setText('Er_C')
        self.tableWidget.setHorizontalHeaderItem(2, item)

        # Eff nodes
        for i_row in range(self.n_Eff):
            # header
            item = QTableWidgetItem()
            item.setText('E{:d}'.format(i_row + 1))
            self.tableWidget.setVerticalHeaderItem(i_row, item)

            # cells
            # item = QCheckBox()
            for i_col in range(self.n_col):
                i_value = self.n_col * i_row + i_col
                item = QDoubleSpinBox()
                item.setValue(values[i_value])
                item.setDecimals(3)
                item.setMinimum(np.finfo(dtype=np.float64).min)
                item.setMaximum(np.finfo(dtype=np.float64).max)
                item.setSingleStep(0.001)
                self.tableWidget.setCellWidget(i_row, i_col, item)

        # nuisance
        for i_row in range(self.n_Eff, self.n_row()):

            nui_ix = list(range(self.n_col * (i_row - self.n_Eff) + 1,
                                np.amin([self.n_nuisance + 1, self.n_col * (i_row - self.n_Eff) + 4])))
            # row header
            item = QTableWidgetItem()
            if len(nui_ix) > 0:
                item.setText('?' + str(nui_ix))
            else:
                item.setText('+/- node')
            self.tableWidget.setVerticalHeaderItem(i_row, item)

            for i_col in range(self.n_col):
                i_value = self.n_col * i_row + i_col
                item = QDoubleSpinBox()
                if i_value < self.n_params():
                    item.setValue(values[i_value])
                else:
                    item.setValue(0.0)
                item.setDecimals(3)
                item.setMinimum(np.finfo(dtype=np.float64).min)
                item.setMaximum(np.finfo(dtype=np.float64).max)
                item.setSingleStep(0.001)
                self.tableWidget.setCellWidget(i_row, i_col, item)

        self.resize_table()
        self.connect_signal()

    def resize_table(self):
        w = self.tableWidget.horizontalHeader().defaultSectionSize()
        h = self.tableWidget.verticalHeader().defaultSectionSize()

        self.tableWidget.setGeometry(QRect(0, 0, w * (self.n_col + 1), h * (self.n_row() + 1)))

    def get_table_size(self):
        print('get_table_size')
        # The method here didn't seem to work
        # w = self.tableWidget.verticalHeader().width() + 4  # +4 seems to be needed
        # for i in range(self.tableWidget.columnCount()):
        #     w += self.tableWidget.columnWidth(i)  # seems to include gridline (on my machine)
        # h = self.tableWidget.horizontalHeader().height() + 4
        # for i in range(self.tableWidget.rowCount()):
        #     h += self.tableWidget.rowHeight(i)

        # add the first column/row one more time
        w = self.tableWidget.horizontalHeader().length() + self.tableWidget.horizontalHeader().sectionSize(0)
        h = self.tableWidget.verticalHeader().length() + self.tableWidget.verticalHeader().sectionSize(0)

        # print(self.tableWidget.horizontalHeader().length())
        # print(self.tableWidget.horizontalHeader().size())
        # print(self.tableWidget.verticalHeader().length())
        # print(self.tableWidget.verticalHeader().size())

        return QSize(w, h)

    def on_value_change(self, value):
        self.value_changed.emit(value)

    def on_stepsize_change(self):
        node = np.zeros(self.n_col, dtype=np.float64)
        for i_col in range(self.n_col):
            item = self.tableWidget.cellWidget(self.tableWidget.rowCount() - 1, i_col)
            node[i_col] = item.value()

        for i_col in range(self.n_col):
            for i_row in range(self.tableWidget.rowCount() - 2):
                i_value = self.n_col * i_row + i_col
                if i_value < self.n_params():
                    item = self.tableWidget.cellWidget(i_row, i_col)
                    item.setSingleStep(node[i_col])

    def connect_signal(self):
        for i_row in range(self.tableWidget.rowCount() - 2):
            for i_col in range(self.tableWidget.columnCount()):
                i_value = self.n_col * i_row + i_col
                if i_value < self.n_params():
                    item = self.tableWidget.cellWidget(i_row, i_col)
                    item.valueChanged.connect(self.on_value_change)
        # step size
        i_row = self.tableWidget.rowCount() - 1
        for i_col in range(self.tableWidget.columnCount()):
            item = self.tableWidget.cellWidget(i_row, i_col)
            item.valueChanged.connect(self.on_stepsize_change)

    def get_value(self):
        theta = []
        for i_row in range(self.tableWidget.rowCount() - 2):
            for i_col in range(self.tableWidget.columnCount()):
                i_value = self.n_col * i_row + i_col
                if i_value < self.n_params():
                    item = self.tableWidget.cellWidget(i_row, i_col)
                    theta.append(item.value())
        # print('get_value')
        print(theta)
        return theta

    def set_value(self, values=[]):
        # print('set_value')
        if len(values) == 0:
            values = [0.0] * self.n_params()

        self.tableWidget.setHorizontalHeaderLabels(['Eff', 'Er_C', 'Er_F'])
        row_header = ['E' + str(i) for i in range(1, self.n_Eff + 1)]
        row_header += ([u"\u03B5" + str(list(range(i, np.amin([self.n_nuisance + 1, i + self.n_col])))) for i in
                        range(1, self.n_nuisance + 1, self.n_col)])
        row_header += ['+/- node', 'Step Size']
        self.tableWidget.setVerticalHeaderLabels(row_header)

        for i_row in range(self.tableWidget.rowCount() - 2):
            for i_col in range(self.tableWidget.columnCount()):
                i_value = self.n_col * i_row + i_col
                if i_value < self.n_params():
                    item = self.tableWidget.cellWidget(i_row, i_col)
                    item.setValue(values[i_value])

        i_row = self.tableWidget.rowCount() - 1
        for i_col in range(self.tableWidget.columnCount()):
            item = self.tableWidget.cellWidget(i_row, i_col)
            item.setValue(0.05)

    def add_node(self):
        # print('add_node')
        # resize
        # add row
        node = np.zeros(self.n_col, dtype=np.float64)
        for i_col in range(self.n_col):
            item = self.tableWidget.cellWidget(self.tableWidget.rowCount() - 2, i_col)
            node[i_col] = item.value()

        old_values = np.array(self.get_value())
        Eff_ix = np.arange(self.n_Eff) * self.n_col
        ix_insert = np.searchsorted(old_values[Eff_ix], node[0], 'left')
        if np.any(old_values[Eff_ix] == node[0]):
            new_values = old_values
            new_values[ix_insert * self.n_col + [1, 2]] = node[[1, 2]]
        else:
            new_values = np.zeros(node.size + old_values.size, dtype=np.float64)
            new_values[:(ix_insert * self.n_col)] = old_values[:(ix_insert * self.n_col)]
            new_values[(ix_insert * self.n_col) + [0, 1, 2]] = node
            new_values[(ix_insert + 1) * self.n_col:] = old_values[(ix_insert * self.n_col):]
            self.n_Eff += 1
            self.add_row(self.tableWidget.rowCount())

        self.set_value(new_values)
        self.resize_table()
        self.connect_signal()
        # print(new_values)
        # self.fill_table(new_values)

    def remove_node(self):
        node = np.zeros(self.n_col, dtype=np.float64)
        for i_col in range(self.n_col):
            item = self.tableWidget.cellWidget(self.tableWidget.rowCount() - 2, i_col)
            node[i_col] = item.value()

        old_values = np.array(self.get_value())
        Eff_ix = np.arange(self.n_Eff) * self.n_col
        ix_remove = np.nonzero(old_values[Eff_ix] == node[0])[0]

        if np.any(ix_remove):
            new_values = old_values
            for i_row in ix_remove[::-1]:
                new_values = np.delete(new_values, np.arange(i_row * self.n_col, (i_row + 1) * self.n_col))
                # print('remove_node')
                # print(old_values)
                # print(new_values)
                # ipdb.set_trace()
                self.remove_row(i_row)
                self.n_Eff -= 1
            self.set_value(new_values)
            self.resize_table()
            self.connect_signal()

    def add_row(self, i_row):
        self.tableWidget.insertRow(i_row)
        for i_col in range(self.n_col):
            item = QDoubleSpinBox()
            item.setDecimals(6)
            item.setMinimum(np.finfo(dtype=np.float64).min)
            item.setMaximum(np.finfo(dtype=np.float64).max)
            item.setSingleStep(0.05)
            item.setKeyboardTracking(False)
            self.tableWidget.setCellWidget(i_row, i_col, item)

    def remove_row(self, i_row):
        self.tableWidget.removeRow(i_row)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
